import pymunk
import pymunk.pygame_util
import pygame
import time
from scipy.spatial import cKDTree
from math import cos, sin, sqrt, atan2
from random import random, uniform, choice
import numpy as np
import matplotlib.pyplot as plt

# Behaviour of panicked crowd
# First version created nov 24 - 2022

### GLOBAL COLOR VARIABLES
backgroundColor  = (46, 52, 64, 255)
wallColor = (76, 86, 106, 255)
particleColor = (216, 222, 233, 255)
deadColor = (191, 97, 106, 255)


space = pymunk.Space()
space._set_sleep_time_threshold(10)
space._set_damping(0.4)



################################################################################





#Calculates angle between two points
def calcAngle(fromPos, toPos):
    deltaY = toPos[1]-fromPos[1]
    deltaX = toPos[0]-fromPos[0]
    return atan2(deltaY,deltaX)

#calculates distance between two points
def getDistance(position,exitposition):
    distance = np.sqrt((position[0]-exitposition[0])**2 + (position[1]-exitposition[1])**2)
    return(distance)




def CalcMeanAngles(neighbours, pointList, particles, nWallPoints, index):
    sinAngles = []
    cosAngles = []
    #Update angles to the mean of all neighbouring balls angles
    for neighbourIndex in neighbours[index]: # neighbours of current boll
        if pointList[neighbourIndex][2]==False:
            sinAngles.append(sin(particles[neighbourIndex-nWallPoints].angle))
            cosAngles.append(cos(particles[neighbourIndex-nWallPoints].angle))
        #meanOfAngles.append(bollar[neighbour].angle) # Append neighbouring balls angle to list
    if len(sinAngles)>0 and len(cosAngles)>0:
        sinAngles = np.mean(sinAngles)
        cosAngles = np.mean(cosAngles)
        meanOfAngles = np.arctan2(sinAngles,cosAngles) + particles[index].noise * np.random.uniform(-1/2,1/2)
    else:
        meanOfAngles = particles[i].angle

    return meanOfAngles

def CalcAvoidAngle(obstacleAvoidance, neighbours, pointList, index, pos):
    xDiff = 0
    yDiff = 0
    if len(neighbours[index]) > 1:
        for obstacleIndex in neighbours[index]:
            if pointList[obstacleIndex][2]==True:
                xDiff+= pos[0]-pointList[obstacleIndex][0]
                yDiff+= pos[1]-pointList[obstacleIndex][1]
            obstacleAngle = atan2((yDiff),(xDiff))
    else:
        obstacleAvoidance = 0.0
        obstacleAngle = np.pi
    if xDiff==0 and yDiff == 0:
        obstacleAvoidance = 0.0
        obstacleAngle = np.pi
    return obstacleAvoidance, obstacleAngle


def CalcSepAngle(separation, neighbours, pointList, index, pos):
    xDiff = 0
    yDiff = 0
    if len(neighbours[index]) > 1:
        for neighbourIndex in neighbours[index]:
            if pointList[neighbourIndex][2]==False and pos[0]!=pointList[neighbourIndex][0] and pos[1]!=pointList[neighbourIndex][1]:
                xDiff+= pos[0]-pointList[neighbourIndex][0]
                yDiff+= pos[1]-pointList[neighbourIndex][1]
        awayAngle = atan2((yDiff),(xDiff))
    else:
        separation = 0.0
        awayAngle = np.pi
    if xDiff==0 and yDiff ==0:
        separation = 0.0
        awayAngle = np.pi
    return separation, awayAngle


######################################################################################



class Wall:
    def __init__(self,startPosition, endPosition, width):
        self.startPosition = startPosition
        self.endPosition = endPosition
        self.length = np.linalg.norm(np.array(self.endPosition) - np.array(self.startPosition)).astype(int)
        self.body = pymunk.Segment(space.static_body, (startPosition), (endPosition), (width))
        self.body.elasticity = 0.0
        self.body.friction = 0.0
        self.body.color = wallColor #Wall color
        space.add(self.body)

    def getBody(self):
        return self.body

    def getPoints(self, resolution):
        nPoints = round(resolution*self.length)
        xPoints = np.linspace(self.startPosition[0],self.endPosition[0],nPoints).astype(list)
        yPoints = np.linspace(self.startPosition[1],self.endPosition[1],nPoints).astype(list)

        isWall = np.full(nPoints,True)

        edgePoints = np.stack((xPoints,yPoints,isWall),axis=1)
        return edgePoints


# A class of the objects we simulate
class Particle:
    def __init__(self,xPos,yPos,angle,setRadius,index,exitPosition):
        self.initialTime = time.time()
        self.body = pymunk.Body(mass=setRadius/np.pi, moment=1)
        self.body.position = xPos, yPos
        self.angle = angle
        self.radius = setRadius
        self.circle = pymunk.Circle(self.body, radius=setRadius)
        self.circle.elasticty = 0.0
        self.circle.friction = 0.0
        self.circle.color = particleColor # Circle color
        self.noise = 0.2
        self.deathConditionIterator = 0
        self.deathProbability = 1/(self.radius**2)
        self.isDead = False
        self.index = index
        self.initialDistance = getDistance([xPos,yPos], exitPosition)
        self.initialPosition = xPos,yPos

        # # Kinetic energy - contribution
        # velo = np.sqrt((self.body._get_velocity()[0])**2 + (self.body._get_velocity()[1])**2)
        # self.initialKE = .5*velo*setRadius/np.pi

        space.add(self.body, self.circle)
    def setAngle(self,newAngle):
            self.angle = newAngle
    def applyForce(self, alignAngle, goalAngle,exitWeight ,sepAngle, separation, obstacleAngle, obstacleAvoidance, time):
        vMax = 50
        force = 100
        a = 0.5 # Alignent constant
        b = separation # Separtion constant
        c = exitWeight*0.1
        d = obstacleAvoidance

        xForce = a*cos(alignAngle)*force + b*cos(sepAngle)*force
        xForce += c*cos(goalAngle)*force + d*cos(obstacleAngle)*force

        yForce = a*sin(alignAngle)*force + b*sin(sepAngle)*force
        yForce += c*sin(goalAngle)*force + d*sin(obstacleAngle)*force


        if  xForce == 0 and yForce ==0:
            pass
        else:
            self.angle = atan2(yForce,xForce)
            self.body.apply_force_at_local_point((xForce, yForce),(0,0))

        # Add vMax restriction
        if abs(self.body._get_velocity()[0]) > vMax or abs(self.body._get_velocity()[1]) > vMax:
            vel = [0,0]
            vel[0] = self.body._get_velocity()[0]
            vel[1] = self.body._get_velocity()[1]
            if vel[0] > vMax:
                vel[0] = vMax
            elif vel[0] <-vMax:
                vel[0] = -vMax
            if vel[1] > vMax:
                vel[1] = vMax
            elif vel[1] <-vMax:
                vel[1] = -vMax
            self.body._set_velocity(vel)

    def checkNeighbours(self):
        touchWall = False
        connectingIndividuals = space.shape_query(self.circle) #
        for individual in connectingIndividuals:
            if individual.shape.body.mass == 0:
                touchWall = True
                break
            if len(connectingIndividuals) > 5:
                self.deathConditionIterator += 1
                if self.deathConditionIterator > 25 and np.random.rand() < self.deathProbability: # if more than 3 consecutive Iterations with more than 5 connecting individuals, kill individual
                    self.isDead = True
                    self.circle.color = deadColor # Circle Dead color
            else:
                self.deathConditionIterator = 0

    def getBody(self):
        return self.body
    def getShape(self):
        return self.circle
##################################################################################################################################
    def Kenergy(self):
        velo = np.sqrt((self.body._get_velocity()[0])**2 + (self.body._get_velocity()[1])**2)
        energy = .5*velo*self.body._get_mass()
        return energy
        self.body._get_

##################################################################################################################################



class App:
    def __init__(self,width,height):
        pygame.init()
        self.screen = pygame.display.set_mode((width, height))
        self.draw_options = pymunk.pygame_util.DrawOptions(self.screen)
        self.running = True
        self.width = width
        self.height = height

        self.exitposition = [width/2,100]
        self.exitDetectionRadius = 350

        #data for plotting dead initial positions
        self.deadPositionListX = []
        self.deadPositionListY = []
        #data for plotting escape time
        self.initialDistances = []
        self.timeSteps = []
        #data for plotting amount of dead, separation and avoidance
        self.amountdead = []
        self.deletedDead = 0
        self.avoidanceList = []
        self.separationList = []


        # This starts the main game loop
    def run(self,particles, obstaclePositions, nrgList):
        self.iteration = 0
        while self.running:
            for event in pygame.event.get():
                if event.type == pygame.QUIT:
                    self.running = False
                    #pygame.image.save(self.screen, 'intro.png')

            pointList = obstaclePositions.copy()
            nWallPoints = len(pointList)

            # Get all positions, we could probably make this faster by making
            # a list of zeros and changing their values instead of appending
            positions = []
            for particle in particles:
                particle.checkNeighbours()
                pos=(particle.getBody()._get_position())
                positions.append([pos[0],pos[1], False])
            pointList.extend(positions)


            #Construct a KDTree and query all neighbours closer than a specified radius
            # The tree contain all points from obstacles and particles, with an extra bool
            # if pointList[i][3] = True the point is from an obstacle, otherwise from a particle
            tree = cKDTree(pointList, boxsize = self.width)
            interactionRadius = 8*particle.radius
            neighbours = tree.query_ball_point(positions,interactionRadius) #
            # List of all parameters for each particle so they are
            # updated synchronmulsity
            foundParameters = False
            updateList = []
            for i in range(len(particles)):
                parameterValues = [0,0,0,0,0,0,0]
                pos = [particles[i].getBody()._get_position()[0],particles[i].getBody()._get_position()[1]]
                if not particles[i].isDead:
                    meanOfAngles = CalcMeanAngles(neighbours, pointList, particles, nWallPoints, i)
                    parameterValues[0] = meanOfAngles

                    # Calculate angle from exitpoint
                    distanceFromExit = getDistance(particles[i].getBody()._get_position(),self.exitposition)
                    if  distanceFromExit < self.exitDetectionRadius: #checks if ball is close to exit and directs ball to exit if so.
                        exitWeight = self.exitDetectionRadius/distanceFromExit       # Creates a normalized value used for weighting the angle update
                        #exitWeight = 0.8
                        exitAngle = calcAngle(particles[i].getBody()._get_position(), [400,-100])
                    else:
                        exitWeight = 0
                        exitAngle = np.pi
                    parameterValues[1] = exitAngle


                    parameterValues[2] = exitWeight
                    # Calculates angle away from other particles
                    separation = 2e-5*len(particles)    # default 2e-5
                    if self.iteration > 100:
                        separation = separation*100/self.iteration
                    separation, awayAngle = CalcSepAngle(separation, neighbours, pointList, i, pos)

                    parameterValues[3] = awayAngle
                    parameterValues[4]= separation


                    # Calculates angle away from obstacles
                    obstacleAvoidance = 0.5
                    if self.iteration > 100:
                        obstacleAvoidance = obstacleAvoidance*100/self.iteration
                    obstacleAvoidance, obstacleAngle = CalcAvoidAngle(obstacleAvoidance, neighbours, pointList, i, pos)
                    parameterValues[5] = obstacleAngle
                    parameterValues[6] = obstacleAvoidance

                    if not foundParameters and separation != 0 and obstacleAvoidance != 0:
                        self.separationList.append(separation)
                        self.avoidanceList.append(obstacleAvoidance)
                        foundParameters = True


                updateList.append(parameterValues)
                    #def applyForce(self, alignAngle, goalAngle,exitWeight ,sepAngle, separation):

            nrgSum = 0
            for iParticle in range(len(particles)):
                if (not particles[iParticle].isDead):
                    particles[iParticle].applyForce(updateList[iParticle][0],updateList[iParticle][1],
                    updateList[iParticle][2],updateList[iParticle][3],updateList[iParticle][4],
                    updateList[iParticle][5],updateList[iParticle][6],self.iteration)
                nrgSum = nrgSum + particles[iParticle].Kenergy()
            nrgList.append(nrgSum)

            #check for dead particles
            deathcounter = 0
            for particle in particles:
                if particle.isDead:
                    deathcounter += 1
            self.amountdead.append(deathcounter+self.deletedDead)


            for particle in particles:
                pos = particle.getBody()._get_position()
                if pos[0] > 850 or pos[0] < 20 or pos[1] > 850 or pos[1] < 20:
                    self.initialDistances.append(particle.initialDistance)
                    self.timeSteps.append(self.iteration)

                    if particle.isDead:
                        self.deadPositionListX.append(particle.initialPosition[0])
                        self.deadPositionListY.append(particle.initialPosition[1])
                        self.deletedDead += 1
                    space.remove(particle.getShape(),particle.getBody())
                    particles.remove(particle)

             # Checking to see if particles have spawn inside the blockade of the central scene and removing them if so - Only applies when central blockade segment is active
            #aX = 350
            #aY = 450
            #for particle in particles:
                #pos = particle.getBody()._get_position()
                #if ((pos[0] > aX) & (pos[0] < aY) & (pos[1] > aX) & (pos[1] < aY)):
                    #space.remove(particle.getShape(),particle.getBody())
                    #particles.remove(particle)

            posit = []
            totalDis = []
            quitDeadNumber = 0
            for i in range(len(particles)):
                posit = particle.getBody()._get_position()
                dis = posit[0], posit[1]
                if particles[i].isDead:
                    quitDeadNumber += 1
            #print(quitDeadNumber,len(particles))
            if quitDeadNumber == len(particles): #if alla remaining particles are dead, quit program.
                break






            if len(particles) <= 0:
                break

            self.screen.fill(backgroundColor)
            space.debug_draw(self.draw_options)
            pygame.display.update()
            #space.step(0.01)
            space.step(0.03)
            time.sleep(0.01)
            self.iteration+=1

        pygame.quit()



def main():
    width = 800
    height = 800
    resolutionOfPoints = .1   #.1 - wall related
    #space.gravity = 0, 100

    window = App(width,height)
    screen = window.screen

    # Create walls around map
    walls = []

    # Central scene
    #aX = 350
    #aY = 450
    #thick = 2
    #walls.append(Wall((aX, aX), (aY, aX), thick))
    #walls.append(Wall((aX, aX), (aX, aY), thick))
    #walls.append(Wall((aY, aX), (aY, aY), thick))
    #walls.append(Wall((aX, aY), (aY, aY), thick))


    #Walls with exit on top Hole with 150
    walls.append(Wall((100, 100), (width/2-50, 100), 4))
    walls.append(Wall((width/2+50, 100), (width-100, 100), 4))
    #Side walls
    walls.append(Wall([100, height-100],[width-100, height-100],4))
    walls.append(Wall([100, 100],[100,height-100],4))
    walls.append(Wall([width-100, 100],[width-100, height-100],4))

    obstaclePositions = []
    for wall in walls:
        obstaclePositions.extend(wall.getPoints(resolutionOfPoints))

    # Create particles with random start position
    # and direction of force vector
    particles = []
    particlesNum = 1200
    for i in range (particlesNum):
        xPos = 580*random()+110
        yPos = 580*random()+110
        currentParticle = Particle(xPos,yPos,uniform(-np.pi,np.pi),np.random.normal(5,0.3),i,window.exitposition)
        particles.append(currentParticle)

     #save particles initial positions:
    initialPositionsX = []
    initialPositionsY = []
    for particle in particles:
        initialPositionsX.append(particle.initialPosition[0])
        initialPositionsY.append(particle.initialPosition[1])


    nrgList = []


    window.run(particles,obstaclePositions, nrgList)



    deadInitialPositionX = window.deadPositionListX #gets the dead particles that were removed
    deadInitialPositionY = window.deadPositionListY
    for particle in particles: #gets the dead particles that remain i nthe window at the end of the run
        if particle.isDead:
            deadInitialPositionX.append(particle.initialPosition[0])
            deadInitialPositionY.append(particle.initialPosition[1])





    fig, ax = plt.subplots()
    plt.scatter(initialPositionsX,initialPositionsY,label='Survived')
    plt.scatter(window.deadPositionListX,window.deadPositionListY,label = 'Died')
    plt.gca().invert_yaxis()
    plt.legend(loc=1)
    plt.title('Initial Positions of Particles')
    #ax.set_facecolor([46/255, 52/255, 64/255,1])
    # plt.show()
    plt.savefig('fig1.png')

    fig, ax = plt.subplots()
    plt.title('Escape Times for Particles based on Starting Position')
    plt.scatter(window.initialDistances,window.timeSteps)
    plt.xlabel('Initial Exit Distance')
    plt.ylabel('Time Steps')
    # plt.show()
    plt.savefig('fig2.png')

    fig, ax = plt.subplots()
    plt.plot(nrgList)
    plt.xlabel('Time Steps')
    plt.ylabel('Kinetic Energy [-]')
    # xx = np.max(nrgList) - .1*np.max(nrgList)
    # yy = len(nrgList) - .1*len(nrgList)
    # s = "Number of particles:" + str(particlesNum)
    # props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    # ax.text(xx, yy, s, transform=ax.transAxes, fontsize=14,
        # verticalalignment='top', bbox=props)
    # ax.text(xx, yy, s, fontdict=None)
    # ax.text(100, 100, "particles")
    # plt.show()
    plt.savefig('fig3.png')
    #print(window.amountdead,range(window.iteration))
    deadplot = [x / particlesNum for x in window.amountdead]
    fig,ax = plt.subplots()
    plt.plot(deadplot)
    plt.title('Accumulative Casualties with respect to Total Population' )
    plt.xlabel('Time Steps')
    plt.ylabel('Casualities [%]')
    plt.savefig('fig4.png')


    fig,ax = plt.subplots()
    plt.title('Change in Panic Parameters over Time')
    plt.plot(window.avoidanceList,label = 'Obstacle Avoidance')
    plt.plot(window.separationList,label = 'Separation')
    plt.xlabel('Time Steps')
    ax.set_yticklabels([])
    plt.legend()
    plt.savefig('fig5.png')


main()
