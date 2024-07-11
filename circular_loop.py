import math
import matplotlib.pyplot as plt

def move_uav_with_known_turn_angle(x, y, V, A, C, time_step, total_time):
    # Initialize lists to store the UAV's path
    path_x = [x]
    path_y = [y]
    
    # Current heading (angle)
    heading = 0  # Assuming initial heading is along the x-axis
    
    for t in range(0, int(total_time / time_step)):
        # Update the heading based on the turn rate and time step
        heading += C * time_step
        
        # Update velocities using the current heading
        x_vel = V * math.cos(heading)
        y_vel = V * math.sin(heading)
        
        # Update positions
        x += x_vel * time_step
        y += y_vel * time_step
        
        # Append new position to path
        path_x.append(x)
        path_y.append(y)
    
    return path_x, path_y

# UAV parameters
initial_x = 0
initial_y = 0
V = 10  # air velocity
A = 2   # turn rate
B = V / A  # turn radius
C = V / B  # turn angle per unit time
time_step = 0.1  # time step for simulation
total_time = 10  # total time for simulation

# Get the UAV's path
path_x, path_y = move_uav_with_known_turn_angle(initial_x, initial_y, V, A, C, time_step, total_time)

# Plot the UAV's path
plt.plot(path_x, path_y, marker='o')
plt.xlabel('X position')
plt.ylabel('Y position')
plt.title('UAV Circular Path with Known Turn Angle')
plt.grid(True)
plt.axis('equal')
plt.show()