import math 
import numpy as np 
from MagneticField import MagneticField
import matplotlib.pyplot as plt 


class GDTS:

    def __init__(self, turn_rate , airspeed):  

        self.turn_radius = airspeed/turn_rate
        self.turn_rate = turn_rate
        self.airspeed = airspeed 
        self.strength = 10 # Strength at the Transmitters Location 
        self.trans_loc = [25,30] # Transmitter Location
        self.total_path_coordinates_x = [] 
        self.total_path_coordinates_y = []
        self.Initialization() 

    def signalstrength(self , pos):
        x , y  = pos[0] , pos[1]
        expont =  - (((x - self.trans_loc[0])**2/0.8) + ((y - self.trans_loc[1])**2/0.8)) 
        curr_strength = self.strength * (math.exp(expont))
        return curr_strength
    
    def magneticField(self , pos_rec , ori_rec , pos_tra , ori_tra):
        obj = MagneticField(pos_rec , ori_rec , pos_tra , ori_tra)
        return obj.magnetic_field 
    

    def Initialization(self):

        self.strengths = [] 

        self.loop = 1
        self.timestep = 0  

        self.del_t = 1/6  
        self.freq = 6

        self.threshold_dist = 0.1

        self.uav_pos = [] 

        x , y = map(float , input("Enter initial X , Y coordinates Seperated by space of Receiver : ").split()) 

        self.initial_rec_loc = [x, y]        
        self.curr_pos = [x , y] 

        self.grad_dir = []
        self.x_vel = 0 
        self.y_vel = 0 

        self.p = [0,0]

        self.w = self.turn_rate 

        self.uav_pos.append(self.curr_pos)

        self.heading = 0 
        self.delta = 2  

        self.IterLoop()  


    def GradientDirection(self , uav_pos , strengths): 
        uav_pos = np.array(uav_pos) 
        A = np.hstack((uav_pos , np.ones((len(strengths) , 1)))) 
        B = np.array(strengths).reshape(-1,1)  
        gradient =  (np.linalg.inv(A.T @ A ) @ A.T ) @ B
        return gradient[:2] 
    
    def IterLoop(self):

        curr_dist = (self.trans_loc[0] - self.uav_pos[-1][0])**2 + (self.trans_loc[1] - self.uav_pos[-1][1])**2  
        self.strengths.append(self.signalstrength(self.curr_pos)) 

        x , y = self.curr_pos

        self.total_path_coordinates_x.append(x)
        self.total_path_coordinates_y.append(y)

        while curr_dist > self.threshold_dist :

            self.timestep += 1   

            self.heading += (self.turn_rate * self.del_t)

            self.x_vel = self.airspeed * np.cos(np.radians(self.heading)) 
            self.y_vel = self.airspeed * np.sin(np.radians(self.heading))  

            x += (self.x_vel *  self.del_t)
            y += (self.y_vel * self.del_t)

            print(f"Current Location of Receiver : {x} , {y} and Target Location is {self.trans_loc[0]} , {self.trans_loc[1]}" )

            self.total_path_coordinates_x.append(x)
            self.total_path_coordinates_y.append(y)

            self.strengths.append(self.signalstrength([x,y])) 
            self.uav_pos.append([x,y]) 

            self.p = np.array([x,y]) - np.array(self.curr_pos) 

            if self.loop == 1 :
                print("I am only Excuting " , self.timestep)  
                if self.timestep > 3 :
                    if self.strengths[self.timestep-3] < self.strengths[self.timestep-2] and self.strengths[self.timestep-2]  == self.strengths[self.timestep-1] and self.strengths[self.timestep-1]  > self.strengths[self.timestep] :
                        self.turn_rate = -self.turn_rate 
                        print("I excuted")
                        self.grad_dir = self.GradientDirection(self.uav_pos , self.strengths) 
                        self.strengths = [] 
                        self.curr_pos = self.uav_pos[-1] 
                        self.uav_pos = [] 
                        self.timestep = 0 
                        self.loop += 1 


                if self.timestep > 2 :
                    if self.strengths[self.timestep-2] < self.strengths[self.timestep-1] and self.strengths[self.timestep-1] > self.strengths[self.timestep] :
                        self.turn_rate = - self.turn_rate 
                        print("I am excuted for 2")
                        self.grad_dir = self.GradientDirection(self.uav_pos , self.strengths)
                        self.strengths = [] 
                        self.curr_pos = self.uav_pos[-1] 
                        self.uav_pos = [] 
                        self.timestep = 0 
                        self.loop += 1 


            else :
                if np.sqrt(np.sum((np.array(self.uav_pos[-1]) - np.array(self.curr_pos))**2)) <= 0.1  :
                    self.turn_rate = - self.turn_rate 
                    self.grad_dir = self.GradientDirection(self.uav_pos , self.strengths) 
                    self.loop += 1 
                    self.timestep = 0
                    self.strengths = [] 
                    self.curr_pos = self.uav_pos[-1] 
                    self.uav_pos = [] 


                if abs(np.degrees(math.acos(np.dot(self.p , self.grad_dir)/(np.linalg.norm(self.p) * np.linalg.norm(self.grad_dir))))) < self.delta : 
                    self.grad_dir = self.GradientDirection(self.uav_pos , self.strengths) 
                    theta = math.atan(self.grad_dir[1]/self.grad_dir[0]) 
                    if (self.heading - theta) * self.turn_rate  >= 0 :
                        self.turn_rate = - self.turn_rate

                    # if np.sqrt(np.sum((np.array(self.uav_pos[-1]) - np.array(self.curr_pos))**2)) <= 0.1  :
                    #     self.turn_rate = - self.turn_rate 
                    
                    self.loop += 1 
                    self.timestep = 0
                    self.strengths = [] 
                    self.curr_pos = self.uav_pos[-1] 
                    self.uav_pos = []      

        self.plotting_path()  


    def plotting_path(self):
        
        plt.figure(figsize=(10,10)) 

        plt.scatter(self.total_path_coordinates_x , self.total_path_coordinates_y , c='r' , s=2)
        plt.scatter([self.trans_loc[0]] , [self.trans_loc[1]] , c='b' , s=3)
        plt.scatter([self.initial_rec_loc[0]] , [self.initial_rec_loc[1]] , c='b' , s=3)
        plt.xlabel('X - coordinates ')
        plt.ylabel('Y - coordinates ')
        plt.title('Total Path Travelled Using GDTS ')
        plt.savefig("GDTS.png")
        plt.show()


if __name__ == '__main__': 
    turn_rate = float(input('Enter the Turn Rate of Receiver: '))
    air_speed = float(input('Enter the Air Speed : '))
    alg = GDTS(turn_rate=turn_rate , airspeed=air_speed) 
