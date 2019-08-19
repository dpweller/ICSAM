from __future__ import division, print_function
from vpython import *

def ShrinkSolidGrowLiquid():#This subroutine changes the sizes of the solid and liquid phases during melting
    solid.pos=solid.pos+vec(0,-(L_box/(2*ceil(m*Hfus/dE))),0) #SizeIncrement1=L_box/(2*ceil(m*Hfus/dE))
    solid.radius=solid.radius-(L_box/(2*ceil(m*Hfus/dE)))
    liquid.pos=liquid.pos+vec(0,0.5*(L_box/(2*ceil(m*Hfus/dE))),0)
    liquid.size=liquid.size+vec(0,(L_box/(2*ceil(m*Hfus/dE))),0)
def ShrinkLiquid():#This subroutine changes the size of the liquid phase during vaporization
    liquid.pos=liquid.pos+vec(0,-0.5*L_box/(2*ceil(m*Hvap/dE)),0) #SizeIncrement2=L_box/(2*ceil(m*Hvap/dE))
    liquid.size=liquid.size+vec(0,-L_box/(2*ceil(m*Hvap/dE)),0)
def MakeNewParticle():#This subroutine creates a new gas particle and adds it to a list of particles
    newparticle = sphere(pos=vector(L_box*(random()-0.5),(L_box-liquid.size.y)*(random()-0.5)+liquid.size.y/2,L_box*(random()-0.5)),radius=0.05, color=color.cyan)
    newparticle.mass = 1
    newparticle.velocity = vector(random()-0.5,random()-0.5,random()-0.5)*2 #the coefficient of 2 is a scaling factor for visual effect
    listOfParticles.append(newparticle)
def ParticleMovementAndCollisions(): #This subroutine moves gas particles and handles particle-wall and particle-particle collisions
    for particle in listOfParticles:
        if Etotal > Eboil: #If all the liquid has evaporated, start increasing the velocity of the gas particles
            particle.velocity=particle.velocity+sqrt(2*dE/particle.mass)*(particle.velocity/mag(particle.velocity))*2E-4 #Increase velocity with increasing energy. The coefficient of 2E-4 is a scaling factor for visual effect
        particle.pos = particle.pos + particle.velocity*0.1 #Update particle position, assume dt=0.1
        if abs(particle.pos.x) >= container.length/2:#Particle-wall collision in x
            particle.velocity.x = - particle.velocity.x
        if abs(particle.pos.y) >= container.height/2 or particle.pos.y <= (liquid.size.y-L_box/2):#Particle-wall collision in y
            particle.velocity.y = - particle.velocity.y
        if abs(particle.pos.z) >= container.width/2:#Particle-wall collision in z
            particle.velocity.z = - particle.velocity.z
    for i in range(0,len(listOfParticles)):#Particle-particle collisions, loop through every particle
        for j in range(i+1,len(listOfParticles)):#loop through every OTHER particle
            diff = listOfParticles[j].pos - listOfParticles[i].pos #displacement vector between two particles
            distance = mag(diff) #magnitude of displacement is distance
            if distance <= listOfParticles[i].radius + listOfParticles[j].radius: #if particles will collide, check their next positions
                nextpos1 = listOfParticles[i].pos + listOfParticles[i].velocity*0.1 #assume dt=0.1
                nextpos2 = listOfParticles[j].pos + listOfParticles[j].velocity*0.1 #assume dt=0.1
                if mag(nextpos2 - nextpos1) < distance: #if they collide with each other, transfer momentum
                    rhat = norm(diff) #unit vector of displacement
                    mv1 = listOfParticles[i].mass*listOfParticles[i].velocity #momentum of first particle
                    mv2 = listOfParticles[j].mass*listOfParticles[j].velocity #momentum of second particle
                    transfer = 2.*dot(listOfParticles[i].mass*mv2-listOfParticles[j].mass*mv1,rhat)/(listOfParticles[i].mass+listOfParticles[j].mass)*rhat #momentum transferred
                    listOfParticles[i].velocity = (mv1 + transfer)/listOfParticles[i].mass
                    listOfParticles[j].velocity = (mv2 - transfer)/listOfParticles[j].mass
