Supplementary Material for Project in the Course FFR120 - Simulation of Complex Systems
----
Emergent Crowd Behaviour Dynamics

Simulating crowd behavior in panic situations using Pygame and Pymunk.

#### Overview

This repository explores the fascinating phenomenon of emergent crowd behavior during emergency scenarios.  It utilizes Pygame for visual rendering and Pymunk to power the realistic physics simulation. The code models a crowd of individuals as they react and attempt to move towards a designated exit under pressure.

##### Key Features

    Individual Agent Behavior: Each crowd member is modeled as an agent with simple rules governing their movement and decision-making.
    Obstacle Interaction: The simulation includes obstacles and walls that agents must navigate around.
    Physics-Based Interactions: Pymunk handles collisions, friction, and other physical forces, making the interactions between agents feel natural.
    Emergent Behaviour: An altered Vicsek model creates crowds through neighbour-aware movement rules for each individual.
##### Prerequisites

    Python 3.x
    Pygame (https://www.pygame.org/)
    Pymunk (https://www.pymunk.org/)
