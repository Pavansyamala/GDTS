import numpy as np 
import pandas as pd
import time 
from datetime import datetime

def setup_data(tran_ori , rec_ori , pr, pt , m_vec = np.array([1, 0, 0])):
    # B_vectors = []

    roll_t = (np.pi/180)*tran_ori[0]  # Roll angle of transmitter 
    pitch_t = (np.pi/180)*tran_ori[1]  # Pitch angle of transmitter 
    yaw_t = (np.pi/180)*tran_ori[2]  # Yaw angle of transmitter 
    
    R_t_to_i = np.array([
        [np.cos(pitch_t)*np.cos(yaw_t), np.sin(roll_t)*np.sin(pitch_t)*np.cos(yaw_t) - np.cos(roll_t)*np.sin(yaw_t), np.cos(roll_t)*np.sin(pitch_t)*np.cos(yaw_t) + np.sin(roll_t)*np.sin(yaw_t)],
        [np.cos(pitch_t)*np.sin(yaw_t), np.sin(roll_t)*np.sin(pitch_t)*np.sin(yaw_t) + np.cos(roll_t)*np.cos(yaw_t), np.cos(roll_t)*np.sin(pitch_t)*np.sin(yaw_t) - np.sin(roll_t)*np.cos(yaw_t)],
        [-np.sin(pitch_t), np.sin(roll_t)*np.cos(pitch_t), np.cos(roll_t)*np.cos(pitch_t)]
    ])

    roll_r = (np.pi/180)*rec_ori[0]  # Roll angle of Reciever 
    pitch_r = (np.pi/180)*rec_ori[1]  # Pitch angle of Reciever
    yaw_r = (np.pi/180)*rec_ori[2]   # Yaw angle of Reciever    
    R_r_to_i = np.array([
        [np.cos(pitch_r)*np.cos(yaw_r), np.sin(roll_r)*np.sin(pitch_r)*np.cos(yaw_r) - np.cos(roll_r)*np.sin(yaw_r), np.cos(roll_r)*np.sin(pitch_r)*np.cos(yaw_r) + np.sin(roll_r)*np.sin(yaw_r)],
        [np.cos(pitch_r)*np.sin(yaw_r), np.sin(roll_r)*np.sin(pitch_r)*np.sin(yaw_r) + np.cos(roll_r)*np.cos(yaw_r), np.cos(roll_r)*np.sin(pitch_r)*np.sin(yaw_r) - np.sin(roll_r)*np.cos(yaw_r)],
        [-np.sin(pitch_r), np.sin(roll_r)*np.cos(pitch_r), np.cos(roll_r)*np.cos(pitch_r)]
    ])

    # B_vectors.append(get_arva_data(pr, pt, R_t_to_i, R_r_to_i, m_vec))
    return  get_arva_data(pr, pt, R_t_to_i, R_r_to_i, m_vec)
    #, R_t_to_i, R_r_to_i

def get_arva_data(pr, pt, R_t_to_i, R_r_to_i, m_vec):
    R_i_to_t = np.transpose(R_t_to_i)
    r = pr - pt
    r = np.dot(R_i_to_t, r)
    A = np.array([2*r[0]**2 - r[1]**2 - r[2]**2 , 3*r[0]*r[1] ,3*r[0]*r[2] ]).reshape(-1,1) 
    Am = np.dot(R_t_to_i, A)
    Am_x = A[0, 0] 
    Am_y = A[1, 0] 
    Am_z = A[2, 0] 
    rd = np.linalg.norm(r)
    H = np.array([(1/(4*np.pi*rd**5))*Am_x, (1/(4*np.pi*rd**5))*Am_y, (1/(4*np.pi*rd**5))*Am_z])
    R_i_to_r = np.transpose(R_r_to_i)
    Hb = np.dot(R_i_to_r, H)
    return Hb 

def orientations():

    roll = [-90 , 90]      
    pitch = [-30 , 30]
    yaw = [-180 , 180] 

    roll_angles = np.linspace(roll[0] , roll[1] , 5)
    pitch_angles = np.linspace(roll[0] , roll[1] , 5)
    yaw_angles = np.linspace(roll[0] , roll[1] , 5) 

    orientations = [] 

    for roll in roll_angles:
        for pitch in pitch_angles:
            for yaw in yaw_angles: 
                orientations.append([roll, pitch, yaw])  
    return orientations 



def transmitterLocations():

    return np.random.randint(low=1 , high= 10 , size=(10,3))


def receiversLocations():
    max_x, min_x = 25, -25
    max_y, min_y = 25, -25
    max_z, min_z = 25, -25

    # Create coordinate grids
    x_vals = np.linspace(min_x, max_x, 21)
    y_vals = np.linspace(min_y, max_y, 21)
    z_vals = np.linspace(min_z, max_z, 21)

    # Generate all combinations of coordinates
    xx, yy, zz = np.meshgrid(x_vals, y_vals, z_vals, indexing='ij')
    locations = np.column_stack((xx.ravel(), yy.ravel(), zz.ravel()))

    return locations


if __name__ == '__main__':

    print("Starting Time of Excution : ", datetime.now().strftime("%Y-%m-%d-%H-%M-%S"))

    start = time.time()
    
    filename = input("Enter the Name with wich you want to save the data generated : ")

    tran_loc = transmitterLocations()
    rec_loc = receiversLocations()  

    reciever_orientations = orientations()
    transmitter_orientations = [[45 , 45 , 45]]

    pos_t , pos_r ,  bstr_vec , ori_rec , ori_tar = [] , [] , [] , [] , []

    count = 1 
    for pt in tran_loc:
        print("Transmitter Location : ", count ) 
        tot_rec_loc = len(rec_loc)
        for pr in rec_loc :
            for i in transmitter_orientations:
                for j in reciever_orientations:

                    b_vec = setup_data(i , j , pr , pt) 
                    pos_t.append(pt)
                    pos_r.append(pr)
                    bstr_vec.append(b_vec)
                    ori_rec.append(j) 
                    ori_tar.append(i)

            print("Rec Loc : " , tot_rec_loc)
            tot_rec_loc -= 1 
        count += 1  


    pos_t = np.array(pos_t)
    pos_r = np.array(pos_r)
    ori_rec = np.array(ori_rec) 
    ori_tar = np.array(ori_tar) 
    bstr_vec = np.array(bstr_vec) 

    pos_trans_wrt_rec = pos_t - pos_r  


    data_dict = {
        "rec_x" : pos_r[:,0] , 
        "rec_y" : pos_r[:,1] ,
        "rec_z" : pos_r[:,2] ,
        "roll_r" : ori_rec[:,0] , 
        "pitch_r" : ori_rec[:,1] ,
        "yaw_r" : ori_rec[:,2] , 
        "mag_x" : bstr_vec[:,0] ,
        "mag_y" : bstr_vec[:,1], 
        "mag_z" : bstr_vec[:,2] ,
        "tar_x" : pos_t[:,0],
        "tar_y" : pos_t[:,1],
        "tar_z" : pos_t[:,2],
        "target_x" : pos_trans_wrt_rec[:,0] , 
        "target_y" : pos_trans_wrt_rec[:,1] , 
        "target_z" : pos_trans_wrt_rec[:,2], 
        "tar_roll" : ori_tar[:,0],
        "tar_pitch" : ori_tar[:,1],
        "tar_yaw" : ori_tar[:,2]
    }

    print(data_dict)

    data = pd.DataFrame(data_dict)
    path_todata = f"D:/DataIITM/{filename}.csv"
    data.to_csv(path_todata , index=False) 

    end = time.time()

    print("Ending Time of Excution : ", datetime.now().strftime("%Y-%m-%d-%H-%M-%S"))

    print("Total Time Taken : " , end-start) 
