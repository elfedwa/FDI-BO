import random
import numpy as np
import csv
import matplotlib.pyplot as plt
from datetime import datetime
from matplotlib.patches import Rectangle
plt.rc('legend', fontsize=10)
def extract_pulses(file_address,sample_size):
       with open(file_address+'profiling.csv', "r") as f:
           numline = (f.readlines())
           print("PRINTING LENGTH OF NUMLINE!!",file_address,len(numline))
           length_of_file=len(numline)
           start=numline[0].split(",")[1]
           finish=numline[length_of_file-1].split(",")[1]
       format = '%H:%M:%S'
       starting_time=datetime.strptime(start.strip(), format)
       day_counter=0
       state_for_incrementing_day_counter=0
       pulses={}
       pulse_counter=0
       for i in range(1,len(numline)-1):
           if file_address == '/Users/awahab/Downloads/NEW_HPC_CAL_DATA_2/SBO-6/':
            print("trubleshoot a ",i, len(pulses), sample_size-1,len(pulses) <= sample_size )
           if len(pulses) <=sample_size:
               if file_address == '/Users/awahab/Downloads/NEW_HPC_CAL_DATA_2/SBO-6/':
                print("trubleshoot b ", i,  len(pulses), sample_size-1,len(pulses) <= sample_size)

               if numline[i].split(",")[0] not in pulses:
                   if (datetime.strptime(numline[i].split(",")[2].strip(), format)-starting_time).total_seconds() >=0 and day_counter==0:
                    pulses[numline[i].split(",")[0]]=[((datetime.strptime(numline[i].split(",")[2].strip(), format)-starting_time).total_seconds())/3600,0]
                    if file_address == '/Users/awahab/Downloads/NEW_HPC_CAL_DATA_2/SBO-6/':
                        print("pulses trouble shoot 1",numline[i].split(",")[0],len(pulses),pulses[numline[i].split(",")[0]])
                   else:
                    a=((datetime.strptime('23:59:59',format) -starting_time).total_seconds())/3600
                    b=((datetime.strptime(numline[i].split(",")[2].strip(), format)-datetime.strptime('0:0:0',format)).total_seconds())/3600
                    pulses[numline[i].split(",")[0]]=[a+b,0]
                    day_counter=1
                    if file_address == '/Users/awahab/Downloads/NEW_HPC_CAL_DATA_2/SBO-6/':

                        print("length pulses trouble shoot 2 with day_counter ",numline[i].split(",")[0],len(pulses),i,pulses[numline[i].split(",")[0]])

               else:
                   #0.0647222
                if file_address == '/Users/awahab/Downloads/NEW_HPC_CAL_DATA_2/SBO-6/':
                    print("identifying id and start and id time",numline[i].split(",")[0],starting_time,datetime.strptime(numline[i].split(",")[2].strip(), format),((datetime.strptime(numline[i].split(",")[2].strip(), format) - starting_time).total_seconds())/3600,day_counter)
                if (datetime.strptime(numline[i].split(",")[2].strip(), format) - starting_time).total_seconds() >0 and day_counter==0:
                 pulses[numline[i].split(",")[0]][1] = ((datetime.strptime(numline[i].split(",")[2].strip(), format) - starting_time).total_seconds())/3600
                 if file_address == '/Users/awahab/Downloads/NEW_HPC_CAL_DATA_2/SBO-6/':

                    print("pulses trouble shoot 3",numline[i].split(",")[0],len(pulses), pulses[numline[i].split(",")[0]][1])

                else:
                 a = ((datetime.strptime('23:59:59', format) - starting_time).total_seconds())/3600
                 b = ((datetime.strptime(numline[i].split(",")[2].strip(), format) - datetime.strptime('0:0:0',
                                                                                                      format)).total_seconds())/3600
                 pulses[numline[i].split(",")[0]][1] = a + b
                 day_counter=1
                 if file_address == '/Users/awahab/Downloads/NEW_HPC_CAL_DATA_2/SBO-6/':

                    print("pulses trouble shoot 4",numline[i].split(",")[0],len(pulses), pulses[numline[i].split(",")[0]])

       print("PULSES LENGTH ADJUSTED",file_address,sample_size,len(pulses))
       if len(pulses)>sample_size:
        pulses.popitem()
        print("POPPED dICTIONARY", len(pulses))
        if file_address == "/Users/awahab/Downloads/NEW_HPC_CAL_DATA_2/new2/SBO-REF-GBR-EI1/":
            print("file name",file_address)
            print("pulses",pulses)
        return [pulses]
       else:
           print("UNPOPPed dICTIONARY", len(pulses))
           if file_address == "/Users/awahab/Downloads/NEW_HPC_CAL_DATA_2/new2/SBO-REF-GBR-EI1/":
               print("file name", file_address)
               print("pulses", pulses)

           return [pulses]


def batch_bin(iter,batch):
 q, mod = divmod(iter, batch)
 if mod > 0:
     bin = q+1
 else:
     bin = q
 return bin

def plotting_profiler(*length):
 fig, ax = plt.subplots(nrows=len(length), sharex=True)
 fig.add_subplot(111, frameon=False)
 finding_feasible_structures={}
 print("TOTAL EXPERIMENTS",len(length))
 reference='S-BO #21'
 plot_reference = 'RS #1'




 cutt_off=-1
 AF_cutt_off=-3
 value_weight=0.8
 time_weight=(1-value_weight)
 new_ratio_with_batch_bin_old = {}
 new_ratio_with_batch_bin_new = {}
 AF={}
 time_calculations_bin_size=30
 best_value_plot={}
 items_that_fall_under_time_bin_criteria=[]

 for iter in range(len(length)):
  finding_feasible_structures[length[iter][0][1]]=len(length[iter][0][2][2])

  new_ratio_with_batch_bin_old[length[iter][0][1]]=((length[iter][0][2][1]['y'] ) / (batch_bin(length[iter][0][2][1]['index'],length[iter][0][3]))) / (length[65][0][2][1]['y'] / batch_bin(length[65][0][2][1]['index'],length[65][0][3]) )  if length[iter][0][2][1]['y'] < cutt_off else 0

  new_ratio_with_batch_bin_new[length[iter][0][1]]=(value_weight*((length[iter][0][2][1]['y'] ) / (length[65][0][2][1]['y'] )) + time_weight*((batch_bin(length[65][0][2][1]['index'],length[65][0][3]) )/(batch_bin(length[iter][0][2][1]['index'],length[iter][0][3]) ))) if length[iter][0][2][1]['y'] < cutt_off else 0

  print("{} FOR THE TRIAL {}, THE MIN WAS FOUND TO BE {} AT ITERATION {} WITH BATCH_SIZE {} SO IT BELONGS TO {} BIN WITH RATIO {}".format(iter,length[iter][0][1],length[iter][0][2][1]['y'],length[iter][0][2][1]['index'],length[iter][0][3],batch_bin(length[iter][0][2][1]['index'],length[iter][0][3]),new_ratio_with_batch_bin_new[length[iter][0][1]]))

 ratio={}
 vec_values={}
 FDI_BO_BEST_VALUE=[]
 FDI_BO_ROP_ARRAY_new=[]
 FDI_BO_ROP_ARRAY_old=[]
 FDI_BO_ROP_ARRAY_AF=[]
 FDI_BO_ROP_ARRAY_time_ROP=[]
 FDI_BO_ROP_ARRAY_time_AF=[]
 FDI_BO_BEST_VALUE_PLOT_samples= {}
 FDI_BO_BEST_VALUE_PLOT_time= {}
 S_BO_BEST_VALUE=[]
 S_BO_ROP_ARRAY_new=[]
 S_BO_ROP_ARRAY_old=[]
 S_BO_ROP_ARRAY_AF = []
 S_BO_ROP_ARRAY_time_ROP = []
 S_BO_ROP_ARRAY_time_AF = []
 S_BO_BEST_VALUE_PLOT_samples= {}
 S_BO_BEST_VALUE_PLOT_time= {}
 TOPK_BO_BEST_VALUE=[]
 TOPK_BO_ROP_ARRAY_new=[]
 TOPK_BO_ROP_ARRAY_old=[]
 TOPK_BO_ROP_ARRAY_AF = []
 TOPK_BO_ROP_ARRAY_time_ROP = []
 TOPK_BO_ROP_ARRAY_time_AF = []
 TOPK_BO_BEST_VALUE_PLOT_samples= {}
 TOPK_BO_BEST_VALUE_PLOT_time= {}
 RANDOM_BEST_VALUE=[]
 RANDOM_ARRAY_new= []
 RANDOM_ARRAY_old= []
 RANDOM_ARRAY_AF = []
 RANDOM_ARRAY_time_ROP = []
 RANDOM_ARRAY_time_AF = []
 RANDOM_BEST_VALUE_PLOT_samples= {}
 RANDOM_BEST_VALUE_PLOT_time={}
 time_campatible_trials = []




 for iter in range(len(length)):
     for i in range(len(length[iter][0][0][0].items())):
         if length[iter][0][1] not in AF:
             AF[length[iter][0][1]] = -1
         print("-AF ERROR BEFOrE", length[iter][0][1], AF[length[iter][0][1]])

         if (length[iter][0][2][0][i]['y']) <=AF_cutt_off:
             AF[length[iter][0][1]]= batch_bin(length[iter][0][2][0][i]['index'],length[iter][0][3])
             print("-AF ERROR AFTER",length[iter][0][1],AF[length[iter][0][1]])

             break
 import copy
 ref_AF=copy.copy(np.float(AF[length[65][0][1]]))

 for iter in range(len(length)):
     AF[length[iter][0][1]]= [ref_AF/np.float(AF[length[iter][0][1]]) if ((AF[length[iter][0][1]]))>0 else 0,AF[length[iter][0][1]]]

 for iter in range(len(length)):
     print("FOR THE TRIAL {} , sample size is {} and batch_size of {} with total BINS {}".format(length[iter][0][1],len(length[iter][0][2][0]),length[iter][0][3],batch_bin(length[iter][0][2][0][len(length[iter][0][2][0])-1]['index'], length[iter][0][3])))
     if batch_bin(length[iter][0][2][0][len(length[iter][0][2][0]) - 1]['index'],
                  length[iter][0][3]) >= time_calculations_bin_size:
         time_campatible_trials.append(length[iter][0][1])
     if 'FDI-BO' in length[iter][0][1]:
         FDI_BO_ROP_ARRAY_new.append(new_ratio_with_batch_bin_new[length[iter][0][1]])
         FDI_BO_ROP_ARRAY_old.append(new_ratio_with_batch_bin_old[length[iter][0][1]])
         FDI_BO_ROP_ARRAY_AF.append(AF[length[iter][0][1]][0])
         FDI_dummy=[]
         for dict_itemss in range(len(length[iter][0][2][3])):
            FDI_dummy.append({'y':length[iter][0][2][3][dict_itemss]['y'],'bin':batch_bin(length[iter][0][2][3][dict_itemss]['index'],
                  length[iter][0][3])})
         FDI_BO_BEST_VALUE_PLOT_samples[length[iter][0][1]]=FDI_dummy
         if batch_bin(length[iter][0][2][0][len(length[iter][0][2][0]) - 1]['index'],
                      length[iter][0][3]) >= time_calculations_bin_size:
             FDI_BO_ROP_ARRAY_time_AF.append(AF[length[iter][0][1]][0])
             FDI_dummy = []
             for dict_itemss in range(len(length[iter][0][2][3])):
                 if batch_bin(length[iter][0][2][3][dict_itemss]['index'],
                              length[iter][0][3]) <= time_calculations_bin_size:
                     FDI_dummy.append({'y': length[iter][0][2][3][dict_itemss]['y'],
                                       'bin': batch_bin(length[iter][0][2][3][dict_itemss]['index'],
                                                        length[iter][0][3])})
                     print("{} with index {} , so bin incremented to {} while the cap is {}".format(length[iter][0][1],length[iter][0][2][3][dict_itemss]['index'],batch_bin(length[iter][0][2][3][dict_itemss]['index'],
                                                                                          length[iter][0][3]),time_calculations_bin_size))

                 else:
                     print("KEY WITH BAD ARRANGEMENT", length[iter][0][1], dict_itemss, length[iter][0][2][3])

                     print("{} terminated due to bin {} so the Best Value Used is {} in bin {}".format(length[iter][0][1],batch_bin(length[iter][0][2][3][dict_itemss]['index'],length[iter][0][3]),FDI_dummy[len(FDI_dummy)-1],FDI_dummy[len(FDI_dummy)-1]['bin']))


                     break
             print("{} the Best Value Used with bin {}".format(length[iter][0][1],
                                                                                                   length[iter][0][3]),
                                                                                               FDI_dummy[
                                                                                                   len(FDI_dummy) - 1],
                                                                                               FDI_dummy[
                                                                                                   len(FDI_dummy) - 1][
                                                                                                   'bin'])

             FDI_BO_ROP_ARRAY_time_ROP.append((value_weight*((FDI_dummy[len(FDI_dummy)-1]['y'] ) / (length[65][0][2][1]['y'] )) + time_weight*((batch_bin(length[65][0][2][1]['index'],length[65][0][3]) )/(FDI_dummy[len(FDI_dummy)-1]['bin'])) if FDI_dummy[len(FDI_dummy)-1]['y'] < cutt_off else 0))

             FDI_BO_BEST_VALUE_PLOT_time[length[iter][0][1]] = FDI_dummy


     elif 'TOPK-BO' in length[iter][0][1]:
         TOPK_BO_ROP_ARRAY_new.append(new_ratio_with_batch_bin_new[length[iter][0][1]])
         TOPK_BO_ROP_ARRAY_old.append(new_ratio_with_batch_bin_old[length[iter][0][1]])
         TOPK_BO_ROP_ARRAY_AF.append(AF[length[iter][0][1]][0])
         TOPK_dummy = []
         for dict_itemss in range(len(length[iter][0][2][3])):
             TOPK_dummy.append({'y': length[iter][0][2][3][dict_itemss]['y'],
                               'bin': batch_bin(length[iter][0][2][3][dict_itemss]['index'],
                                                length[iter][0][3])})
         TOPK_BO_BEST_VALUE_PLOT_samples[length[iter][0][1]] = TOPK_dummy
         if batch_bin(length[iter][0][2][0][len(length[iter][0][2][0]) - 1]['index'],
                      length[iter][0][3]) >= time_calculations_bin_size:
             # TOPK_BO_ROP_ARRAY_time_ROP.append(new_ratio_with_batch_bin_new[length[iter][0][1]])
             print("WHY AF IS NEGATIVE",AF[length[iter][0][1]])
             TOPK_BO_ROP_ARRAY_time_AF.append(AF[length[iter][0][1]][0])
             TOPK_dummy = []
             for dict_itemss in range(len(length[iter][0][2][3])):
                 if batch_bin(length[iter][0][2][3][dict_itemss]['index'],
                              length[iter][0][3]) <= time_calculations_bin_size:
                     TOPK_dummy.append({'y': length[iter][0][2][3][dict_itemss]['y'],
                                       'bin': batch_bin(length[iter][0][2][3][dict_itemss]['index'],
                                                        length[iter][0][3])})
                 else:
                     print("{} terminated due to {} bin".format(length[iter][0][1],
                                                                batch_bin(length[iter][0][2][3][dict_itemss]['index'],
                                                                          length[iter][0][3])))
                     break
             print("{} the Best Value Used with bin {}".format(length[iter][0][1],
                                                               length[iter][0][3]),
                   TOPK_dummy[
                       len(TOPK_dummy) - 1],
                   TOPK_dummy[
                       len(TOPK_dummy) - 1][
                       'bin'])
             TOPK_BO_ROP_ARRAY_time_ROP.append((value_weight*((TOPK_dummy[len(TOPK_dummy)-1]['y'] ) / (length[65][0][2][1]['y'] )) + time_weight*((batch_bin(length[65][0][2][1]['index'],length[65][0][3]) )/(TOPK_dummy[len(TOPK_dummy)-1]['bin'])) if TOPK_dummy[len(TOPK_dummy)-1]['y'] < cutt_off else 0))

             TOPK_BO_BEST_VALUE_PLOT_time[length[iter][0][1]] = TOPK_dummy


     elif 'S-BO' in length[iter][0][1]:
         S_BO_ROP_ARRAY_new.append(new_ratio_with_batch_bin_new[length[iter][0][1]])
         S_BO_ROP_ARRAY_old.append(new_ratio_with_batch_bin_old[length[iter][0][1]])
         S_BO_ROP_ARRAY_AF.append(AF[length[iter][0][1]][0])
         SBO_dummy = []
         for dict_itemss in range(len(length[iter][0][2][3])):
             SBO_dummy.append({'y': length[iter][0][2][3][dict_itemss]['y'],
                                'bin': batch_bin(length[iter][0][2][3][dict_itemss]['index'],
                                                 length[iter][0][3])})
         S_BO_BEST_VALUE_PLOT_samples[length[iter][0][1]] = SBO_dummy
         if batch_bin(length[iter][0][2][0][len(length[iter][0][2][0]) - 1]['index'],
                      length[iter][0][3]) >= time_calculations_bin_size:
             S_BO_ROP_ARRAY_time_AF.append(AF[length[iter][0][1]][0])
             SBO_dummy = []
             for dict_itemss in range(len(length[iter][0][2][3])):
                 if batch_bin(length[iter][0][2][3][dict_itemss]['index'],
                              length[iter][0][3]) <= time_calculations_bin_size:
                     SBO_dummy.append({'y': length[iter][0][2][3][dict_itemss]['y'],
                                        'bin': batch_bin(length[iter][0][2][3][dict_itemss]['index'],
                                                         length[iter][0][3])})
                 else:
                     print("{} terminated due to {} bin".format(length[iter][0][1],
                                                                batch_bin(length[iter][0][2][3][dict_itemss]['index'],
                                                                          length[iter][0][3])))
                     break
             print("{} the Best Value Used with bin {}".format(length[iter][0][1],
                                                               length[iter][0][3]),
                   SBO_dummy[
                       len(SBO_dummy) - 1],
                   SBO_dummy[
                       len(SBO_dummy) - 1][
                       'bin'])
             S_BO_ROP_ARRAY_time_ROP.append((value_weight * (
                         (SBO_dummy[len(SBO_dummy) - 1]['y']) / (length[65][0][2][1]['y'])) + time_weight * ((
                                                                                                                   batch_bin(
                                                                                                                       length[
                                                                                                                           65][
                                                                                                                           0][
                                                                                                                           2][
                                                                                                                           1][
                                                                                                                           'index'],
                                                                                                                       length[
                                                                                                                           65][
                                                                                                                           0][
                                                                                                                           3])) / (
                                                                                                               SBO_dummy[
                                                                                                                   len(SBO_dummy) - 1][
                                                                                                                   'bin'])) if
                                                SBO_dummy[len(SBO_dummy) - 1]['y'] < cutt_off else 0))

             S_BO_BEST_VALUE_PLOT_time[length[iter][0][1]] = SBO_dummy



     elif 'RS' in length[iter][0][1]:
         RANDOM_ARRAY_new.append(new_ratio_with_batch_bin_new[length[iter][0][1]])
         RANDOM_ARRAY_old.append(new_ratio_with_batch_bin_old[length[iter][0][1]])
         RANDOM_ARRAY_AF.append(AF[length[iter][0][1]][0])
         RANDOM_dummy = []
         for dict_itemss in range(len(length[iter][0][2][3])):
             RANDOM_dummy.append({'y': length[iter][0][2][3][dict_itemss]['y'],
                               'bin': batch_bin(length[iter][0][2][3][dict_itemss]['index'],
                                                length[iter][0][3])})
         RANDOM_BEST_VALUE_PLOT_samples[length[iter][0][1]] = RANDOM_dummy
         if batch_bin(length[iter][0][2][0][len(length[iter][0][2][0]) - 1]['index'],
                      length[iter][0][3]) >= time_calculations_bin_size:
             RANDOM_ARRAY_time_AF.append(AF[length[iter][0][1]][0])
             RANDOM_dummy = []
             for dict_itemss in range(len(length[iter][0][2][3])):
                 if batch_bin(length[iter][0][2][3][dict_itemss]['index'],
                              length[iter][0][3]) <= time_calculations_bin_size:
                     RANDOM_dummy.append({'y': length[iter][0][2][3][dict_itemss]['y'],
                                        'bin': batch_bin(length[iter][0][2][3][dict_itemss]['index'],
                                                         length[iter][0][3])})
                 else:
                     print("{} terminated due to {} bin".format(length[iter][0][1],
                                                                batch_bin(length[iter][0][2][3][dict_itemss]['index'],
                                                                          length[iter][0][3])))
                     break
             print("{} the Best Value Used with bin {}".format(length[iter][0][1],
                                                               length[iter][0][3]),
                   RANDOM_dummy[
                       len(RANDOM_dummy) - 1],
                   RANDOM_dummy[
                       len(RANDOM_dummy) - 1][
                       'bin'])
             RANDOM_ARRAY_time_ROP.append((value_weight * (
                     (RANDOM_dummy[len(RANDOM_dummy) - 1]['y']) / (length[65][0][2][1]['y'])) + time_weight * ((
                                                                                                             batch_bin(
                                                                                                                 length[
                                                                                                                     65][
                                                                                                                     0][
                                                                                                                     2][
                                                                                                                     1][
                                                                                                                     'index'],
                                                                                                                 length[
                                                                                                                     65][
                                                                                                                     0][
                                                                                                                     3])) / (
                                                                                                             RANDOM_dummy[
                                                                                                                 len(RANDOM_dummy) - 1][
                                                                                                                 'bin'])) if
                                             RANDOM_dummy[len(RANDOM_dummy) - 1]['y'] < cutt_off else 0))
             RANDOM_BEST_VALUE_PLOT_time[length[iter][0][1]] = RANDOM_dummy
     vecx=[]
     vecy=[]
     for j in range(len(length[iter][0][2][3])):
         # print("j is ", length[iter][0][2][3][j])
         x = list(length[iter][0][0][0].values())[length[iter][0][2][3][j]['index']-1][1]

         vecx.append(list(length[iter][0][0][0].values())[length[iter][0][2][3][j]['index']-1][1])
         vecy.append(length[iter][0][2][3][j]['y'])

     if (vecy[len(vecy) - 1]) > -1:
         ratio[length[iter][0][1]] = 0
     else:
         ratio[length[iter][0][1]] = abs(round(((vecy[len(vecy) - 1]) / vecx[len(vecx) - 1]) * 100, 2))

     vec_values[length[iter][0][1]] = {'x': vecx, 'y': vecy}

 FDI_only_BEST=[]
 TOPK_only_BEST=[]
 SBO_only_BEST=[]
 RS_only_BEST=[]

 for key in FDI_BO_BEST_VALUE_PLOT_samples.keys():
     # print("FDI key {} best value {}".format(key,list(FDI_BO_BEST_VALUE_PLOT_samples[key])[len(FDI_BO_BEST_VALUE_PLOT_samples[key])-1]['y']))
     FDI_only_BEST.append(list(FDI_BO_BEST_VALUE_PLOT_samples[key])[len(FDI_BO_BEST_VALUE_PLOT_samples[key])-1]['y'])
 for key in TOPK_BO_BEST_VALUE_PLOT_samples.keys():
     # print("FDI key {} best value {}".format(key,list(FDI_BO_BEST_VALUE_PLOT_samples[key])[len(FDI_BO_BEST_VALUE_PLOT_samples[key])-1]['y']))
     TOPK_only_BEST.append(list(TOPK_BO_BEST_VALUE_PLOT_samples[key])[len(TOPK_BO_BEST_VALUE_PLOT_samples[key])-1]['y'])
 for key in S_BO_BEST_VALUE_PLOT_samples.keys():
     # print("FDI key {} best value {}".format(key,list(FDI_BO_BEST_VALUE_PLOT_samples[key])[len(FDI_BO_BEST_VALUE_PLOT_samples[key])-1]['y']))
     SBO_only_BEST.append(list(S_BO_BEST_VALUE_PLOT_samples[key])[len(S_BO_BEST_VALUE_PLOT_samples[key])-1]['y'])
 for key in RANDOM_BEST_VALUE_PLOT_samples.keys():
     # print("FDI key {} best value {}".format(key,list(FDI_BO_BEST_VALUE_PLOT_samples[key])[len(FDI_BO_BEST_VALUE_PLOT_samples[key])-1]['y']))
     RS_only_BEST.append(list(RANDOM_BEST_VALUE_PLOT_samples[key])[len(RANDOM_BEST_VALUE_PLOT_samples[key])-1]['y'])

 if sample_size < 200:
     print("AVERAGE BEST VAL FOR FDI-BO", np.round(np.mean(FDI_only_BEST), 3))
     print("AVERAGE BEST VAL FOR TOPK-BO", np.round(np.mean(TOPK_only_BEST), 3))
     print("AVERAGE BEST VAL FOR S-BO", np.round(np.mean(SBO_only_BEST), 3))
     print("AVERAGE BEST VAL FOR RANDOM", np.round(np.mean(RS_only_BEST), 3))
     print('###########################################################################################')
     print("AVERAGE new ROP FOR TOPK-BO",np.round(np.mean(TOPK_BO_ROP_ARRAY_new),3),value_weight,np.round(time_weight,2))
     print("AVERAGE new ROP FOR FDI-BO",np.round(np.mean(FDI_BO_ROP_ARRAY_new),3),value_weight,np.round(time_weight,2))
     print("AVERAGE new ROP FOR S-BO",np.round(np.mean(S_BO_ROP_ARRAY_new),3),value_weight,np.round(time_weight,2))
     print("AVERAGE new ROP FOR RANDOM",np.round(np.mean(RANDOM_ARRAY_new),3),value_weight,np.round(time_weight,2))

     print('###########################################################################################')
     print("AVERAGE AF FOR TOPK-BO",np.round(np.mean(TOPK_BO_ROP_ARRAY_AF),3))
     print("AVERAGE AF FOR FDI-BO",np.round(np.mean(FDI_BO_ROP_ARRAY_AF),3))
     print("AVERAGE AF FOR S-BO",np.round(np.mean(S_BO_ROP_ARRAY_AF),3))
     print("AVERAGE AF FOR RANDOM",np.round(np.mean(RANDOM_ARRAY_AF),3))
     print('###########################################################################################')

 if sample_size > 200:
     print("AVERAGE BEST VAL FOR FDI-BO", np.round(np.mean(FDI_only_BEST), 3))
     print("AVERAGE BEST VAL FOR TOPK-BO", np.round(np.mean(TOPK_only_BEST), 3))
     print("AVERAGE BEST VAL FOR S-BO", np.round(np.mean(SBO_only_BEST), 3))
     print("AVERAGE BEST VAL FOR RANDOM", np.round(np.mean(RS_only_BEST), 3))
     print('###########################################################################################')
     print("TOTAL STRUCTURES & TIME COMPATIBLE TRIALS ARE",len(time_campatible_trials),time_campatible_trials)
     print('###########################################################################################')
     print("AVERAGE TIME ROP FOR TOPK-BO", len(TOPK_BO_ROP_ARRAY_time_ROP),np.round(np.mean(TOPK_BO_ROP_ARRAY_time_ROP), 3))
     print("AVERAGE TIME ROP FOR FDI-BO", len(FDI_BO_ROP_ARRAY_time_ROP),np.round(np.mean(FDI_BO_ROP_ARRAY_time_ROP), 3))
     print("AVERAGE TIME ROP FOR S-BO", len(S_BO_ROP_ARRAY_time_ROP),np.round(np.mean(S_BO_ROP_ARRAY_time_ROP), 3))
     print("AVERAGE TIME ROP FOR RANDOM", len(RANDOM_ARRAY_time_ROP),np.round(np.mean(RANDOM_ARRAY_time_ROP), 3))
     print('###########################################################################################')
     print("AVERAGE TIME AF FOR TOPK-BO", len(TOPK_BO_ROP_ARRAY_time_AF),np.round(np.mean(TOPK_BO_ROP_ARRAY_time_AF), 3))
     print("AVERAGE TIME AF FOR FDI-BO", len(FDI_BO_ROP_ARRAY_time_AF),np.round(np.mean(FDI_BO_ROP_ARRAY_time_AF), 3))
     print("AVERAGE TIME AF FOR S-BO", len(S_BO_ROP_ARRAY_time_AF),np.round(np.mean(S_BO_ROP_ARRAY_time_AF), 3))
     print("AVERAGE TIME AF FOR RANDOM", len(RANDOM_ARRAY_time_AF),np.round(np.mean(RANDOM_ARRAY_time_AF), 3))
 import copy
 reference_copy=copy.copy(ratio[reference])
 print("ratio of reference",reference_copy)


#################################################################################################### PLOTTING ROP AVERAGE ###############################################################################################################################################################
 if sample_size < 200:
     plt.figure()


     plt.bar(x='FDI-BO',

             height=round(np.round(np.mean(FDI_BO_ROP_ARRAY_new),3), 2), alpha=0.6, color='darkgray', hatch='//', capsize=12)

     plt.bar(x='S-BO',

             height=round(np.round(np.mean(S_BO_ROP_ARRAY_new),3), 2), color='gray', alpha=0.6, hatch='xx', capsize=12)

     plt.bar(x='TOPK-BO',
             height=round(np.round(np.mean(TOPK_BO_ROP_ARRAY_new),3), 2), color='dimgray', alpha=0.6, hatch='..', capsize=12)

     plt.bar(x='Random',
             height=round(np.round(np.mean(RANDOM_ARRAY_new),3), 2), color='dimgray', alpha=0.6, hatch='./.', capsize=12)

     plt.xticks(rotation=20)
     plt.ylabel('Average Relative Overall Performance')
     plt.title('α={} and β ={} for Fixed Samples Trials'.format(value_weight,np.round(time_weight,2)))
     # plt.savefig('/Users/awahab/Documents/QEERI/IMMI_REVISION/IMMI_NEW_EXPERIMENTS_FIGURES/ROP_Fixed_Samples.png')
     plt.tight_layout()
 ######################################################################################################################################## ###### ACCELERATION FACTOR ######
##################################################samp
 if sample_size<200:
     plt.figure()
     ###### ACCELERATION FACTOR ######



     import statistics

     plt.bar(x='FDI-BO',

             height=round(np.round(np.mean(FDI_BO_ROP_ARRAY_AF),3), 2), alpha=0.6, color='darkgray', hatch='//', capsize=12)
     plt.bar(x='S-BO',
             height=round(np.round(np.mean(S_BO_ROP_ARRAY_AF), 3), 2), color='gray', alpha=0.6, hatch='xx', capsize=12
             )


     plt.bar(x='TOPK-BO',
             height=round(np.round(np.mean(TOPK_BO_ROP_ARRAY_AF),3), 2), color='dimgray', alpha=0.6, hatch='..', capsize=12)

     plt.bar(x='RS',
             height=round(np.round(np.mean(RANDOM_ARRAY_AF), 3), 2), color='dimgray', alpha=0.6, hatch='./.', capsize=12
             )

     plt.xticks(rotation=20)
     plt.ylabel('Average Acceleration Factor')
     plt.title('Fixed Samples Trials with {} Cut-off'.format(AF_cutt_off))
     plt.savefig('/Users/awahab/Documents/QEERI/IMMI_REVISION/IMMI_NEW_EXPERIMENTS_FIGURES/AF_Corrected_for_lower_bound.png')

 ############################################################################################################################# TIME ROP######################################################################################################################################
 if sample_size>200:
     plt.figure()


     plt.bar(x='FDI-BO',

             height=round(np.round(np.mean(FDI_BO_ROP_ARRAY_time_ROP),3), 2), alpha=0.6, color='darkgray', hatch='//', capsize=12)

     plt.bar(x='S-BO',

             height=round(np.round(np.mean(S_BO_ROP_ARRAY_time_ROP),3), 2), color='gray', alpha=0.6, hatch='xx', capsize=12)

     plt.bar(x='TOPK-BO',
             height=round(np.round(np.mean(TOPK_BO_ROP_ARRAY_time_ROP),3), 2), color='dimgray', alpha=0.6, hatch='..', capsize=12)

     plt.bar(x='RS',
             height=round(np.round(np.mean(RANDOM_ARRAY_time_ROP),3), 2), color='dimgray', alpha=0.6, hatch='./.', capsize=12)

     plt.xticks(rotation=20)
     plt.ylabel('Average Relative Overall Performance')
     plt.title('α={} and β ={} for Fixed Time Trials'.format(value_weight,np.round(time_weight,2)))
     # plt.savefig('/Users/awahab/Documents/QEERI/IMMI_REVISION/IMMI_NEW_EXPERIMENTS_FIGURES/ROP_Fixed_Time.png')

     plt.tight_layout()
 ############################################################################################################################# TIME AF ######################################################################################################################################
 if sample_size > 200:
     plt.figure()

     plt.bar(x='FDI-BO',

             height=round(np.round(np.mean(FDI_BO_ROP_ARRAY_time_AF), 3), 2), alpha=0.6, color='darkgray', hatch='//',
             capsize=12)
     plt.bar(x='S-BO',
             height=round(np.round(np.mean(S_BO_ROP_ARRAY_time_AF), 3), 2), color='gray', alpha=0.6, hatch='xx',
             capsize=12
             )


     plt.bar(x='TOPK-BO',
             height=round(np.round(np.mean(TOPK_BO_ROP_ARRAY_time_AF), 3), 2), color='dimgray', alpha=0.6, hatch='..',
             capsize=12)

     plt.bar(x='RS',

             height=round(np.round(np.mean(RANDOM_ARRAY_time_AF), 3), 2), color='dimgray', alpha=0.6, hatch='./.',
             capsize=12)

     plt.xticks(rotation=20)
     plt.ylabel('Average Acceleration Factor')
     plt.title('Fixed Time Trials with {} Cut-off'.format(AF_cutt_off))
     plt.savefig('/Users/awahab/Documents/QEERI/IMMI_REVISION/IMMI_NEW_EXPERIMENTS_FIGURES/AF_Fixed_Time.png')

     plt.tight_layout()
 ################################################################################################################### PLOTTONG BEST_VALUES REACHED ################################################################################################################################################
 plt.figure()
 fd=0
 for ii in range(len(FDI_BO_BEST_VALUE_PLOT_samples)):
     for key in FDI_BO_BEST_VALUE_PLOT_samples.keys():
         y=[list(FDI_BO_BEST_VALUE_PLOT_samples[key])[values]['y'] for values in range(len(FDI_BO_BEST_VALUE_PLOT_samples[key]))]
         x=[list(FDI_BO_BEST_VALUE_PLOT_samples[key])[values]['bin'] for values in range(len(FDI_BO_BEST_VALUE_PLOT_samples[key]))]

         plt.scatter(x,y,color='r',s=0.5)
         if fd == 0:
             fd=1
             plt.plot(x, y, linewidth='0.5',color='r',label='FDI-BO')
         else:
             plt.plot(x, y, linewidth='0.5',color='r')
         for hmm in range(len(y)):
             if y[hmm] <= -4:
                 plt.annotate((str(round(y[hmm], 3)) + ': ' + key), xy=(
                     x[hmm], y[hmm]), xytext=(x[hmm] + 1, y[hmm] - 0.5), size=18, ha='center', \
                              )
 tk=0
 for ij in range(len(TOPK_BO_BEST_VALUE_PLOT_samples)):
     for key in TOPK_BO_BEST_VALUE_PLOT_samples.keys():
         y=[list(TOPK_BO_BEST_VALUE_PLOT_samples[key])[values]['y'] for values in range(len(TOPK_BO_BEST_VALUE_PLOT_samples[key]))]
         x=[list(TOPK_BO_BEST_VALUE_PLOT_samples[key])[values]['bin'] for values in range(len(TOPK_BO_BEST_VALUE_PLOT_samples[key]))]
         plt.scatter(x,y,color='blue',s=0.5)
         if tk ==0:
             plt.plot(x, y,color='blue', linewidth='0.5',label='TOPK-BO')
             tk=1
         else:
             plt.plot(x, y, color='blue', linewidth='0.5')
         for hmm in range(len(y)):
             if y[hmm] <= -5:
                 plt.annotate((str(round(y[hmm], 3)) + ': ' + key), xy=(
                     x[hmm], y[hmm]), xytext=(x[hmm] + 1, y[hmm] - 0.5), size=18, ha='center', \
                              )

 sb=0
 for ik in range(len(S_BO_BEST_VALUE_PLOT_samples)):
     for key in S_BO_BEST_VALUE_PLOT_samples.keys():
         y=[list(S_BO_BEST_VALUE_PLOT_samples[key])[values]['y'] for values in range(len(S_BO_BEST_VALUE_PLOT_samples[key]))]
         x=[list(S_BO_BEST_VALUE_PLOT_samples[key])[values]['bin'] for values in range(len(S_BO_BEST_VALUE_PLOT_samples[key]))]
         plt.scatter(x,y,color='green',s=0.5)
         if sb == 0:
            plt.plot(x, y,color='green', linewidth='0.5',label='S-BO')
            sb=1
         else:
            plt.plot(x, y, color='green', linewidth='0.5')
         for hmm in range(len(y)):
             if y[hmm] <= -5:
                 plt.annotate((str(round(y[hmm], 3)) + ': ' + key), xy=(
                     x[hmm], y[hmm]), xytext=(x[hmm] + 1, y[hmm] - 0.5), size=18, ha='center', \
                              )
 rnd=0
 for il in range(len(RANDOM_BEST_VALUE_PLOT_samples)):
     for key in RANDOM_BEST_VALUE_PLOT_samples.keys():
         y=[list(RANDOM_BEST_VALUE_PLOT_samples[key])[values]['y'] for values in range(len(RANDOM_BEST_VALUE_PLOT_samples[key]))]
         x=[list(RANDOM_BEST_VALUE_PLOT_samples[key])[values]['bin'] for values in range(len(RANDOM_BEST_VALUE_PLOT_samples[key]))]
         plt.scatter(x,y,color='brown',s=0.5)
         if rnd ==0:
            rnd=1
            plt.plot(x, y,color='black',label='RS',linewidth='0.5')
         else:
            plt.plot(x, y, color='black',linewidth='0.5')
         for hmm in range(len(y)):
             if y[hmm] <= -5:
                 plt.annotate((str(round(y[hmm], 3)) + ': ' + key), xy=(
                     x[hmm], y[hmm]), xytext=(x[hmm] + 1, y[hmm] - 0.5), size=18, ha='center', \
                              )

 plt.legend(loc='upper right')
 plt.xticks(fontsize=15)
 plt.yticks(fontsize=15)

 plt.ylabel('ΔHmix (meV/ion)',fontsize = 20)
 plt.xlabel('Batch Cycle',fontsize = 20)

 ################################################################################################################### PLOTTONG PROFILES ################################################################################################################################################

 track_which_ones=[]
 for k in range(5):
     fig, ax = plt.subplots(nrows=20, sharex=True, figsize=(7, 9))
     fig.add_subplot(111, frameon=False)

     # plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.5)
     print("length of length", len(length))
     plt.ylabel('ΔHmix (meV/ion)', fontsize=10, labelpad=15)
     plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
     universal_counter = 0
     FDI_budget = 5
     FDI_label=0
     TOKK_budget = 5
     TOKK_label=0
     SBO_Budget = 5
     SBO_label=0
     RS_Budget = 5
     RS_label=0
     for iter in range(len(length)):
         if 'FDI-BO' in length[iter][0][1] and FDI_budget>=1 and (length[iter][0][1] not in track_which_ones):
            for i, (key, value) in enumerate(length[iter][0][0][0].items()):
                if value[1] is not 0:
                    if FDI_label==0:
                         ax[universal_counter].add_patch(
                             Rectangle((value[0], 0), value[1] - value[0], length[iter][0][2][0][i]['y'], lw=1, fill=None, color='green',label="FDI-BO #" + str(1+(k*5)) + " - #" + str(5+(k*5))))
                         ax[universal_counter].legend(loc='best', markerscale=0.25, handlelength=0.25)
                         FDI_label=1
                    else:
                        ax[universal_counter].add_patch(
                            Rectangle((value[0], 0), value[1] - value[0], length[iter][0][2][0][i]['y'], lw=1,
                                      fill=None, color='green'))
            x = vec_values[length[iter][0][1]]['x'][len(vec_values[length[iter][0][1]]['x']) - 1]
            ax[universal_counter].axvline(x, color='cyan', linestyle='--')
            ax[universal_counter].set_ylim([-8,8])
            universal_counter+=1
            track_which_ones.append(length[iter][0][1])
            FDI_budget-=1
         print("universal_counter FDI-BO",universal_counter)
     for iter in range(len(length)):

         if 'TOPK-BO' in length[iter][0][1] and TOKK_budget>=1  and (length[iter][0][1] not in track_which_ones):
            for i, (key, value) in enumerate(length[iter][0][0][0].items()):
                if value[1] is not 0:
                    if TOKK_label==0:
                        ax[universal_counter].add_patch(
                            Rectangle((value[0], 0), value[1] - value[0], length[iter][0][2][0][i]['y'], lw=1,
                                      fill=None, color='black',label="TOPK-BO #" + str(1+(k*5)) + " - #" + str(5+(k*5))))
                        ax[universal_counter].legend(loc='best', markerscale=0.25, handlelength=0.25)
                        TOKK_label=1

                    else:
                         ax[universal_counter].add_patch(
                             Rectangle((value[0], 0), value[1] - value[0], length[iter][0][2][0][i]['y'], lw=1, fill=None, color='black'))
            x = vec_values[length[iter][0][1]]['x'][len(vec_values[length[iter][0][1]]['x']) - 1]
            ax[universal_counter].axvline(x, color='cyan', linestyle='--')
            ax[universal_counter].set_ylim([-10, 10])
            universal_counter += 1
            track_which_ones.append(length[iter][0][1])
            TOKK_budget -= 1
     for iter in range(len(length)):

         if 'S-BO' in length[iter][0][1] and SBO_Budget>=1  and (length[iter][0][1] not in track_which_ones):
            for i, (key, value) in enumerate(length[iter][0][0][0].items()):
                if value[1] is not 0:
                    if SBO_label == 0:
                         ax[universal_counter].add_patch(
                             Rectangle((value[0], 0), value[1] - value[0], length[iter][0][2][0][i]['y'], lw=1, fill=None, color='blue',label="S-BO #" + str(1+(k*5)) + " - #" + str(5+(k*5))))
                         ax[universal_counter].legend(loc='upper right', markerscale=0.25, handlelength=0.25)
                         SBO_label = 1

                    else:
                        ax[universal_counter].add_patch(
                            Rectangle((value[0], 0), value[1] - value[0], length[iter][0][2][0][i]['y'], lw=1,
                                      fill=None, color='blue'))
            x = vec_values[length[iter][0][1]]['x'][len(vec_values[length[iter][0][1]]['x']) - 1]
            ax[universal_counter].axvline(x, color='cyan', linestyle='--')
            ax[universal_counter].set_ylim([-10, 10])
            universal_counter += 1
            track_which_ones.append(length[iter][0][1])
            SBO_Budget -= 1
     for iter in range(len(length)):

         if 'RS' in length[iter][0][1] and RS_Budget>=1  and (length[iter][0][1] not in track_which_ones):
            for i, (key, value) in enumerate(length[iter][0][0][0].items()):
                if value[1] is not 0:
                    if RS_label == 0:
                         ax[universal_counter].add_patch(
                             Rectangle((value[0], 0), value[1] - value[0], length[iter][0][2][0][i]['y'], lw=1, fill=None, color='red',label="RS #" + str(1+(k*5)) + " - #" + str(5+(k*5))))
                         ax[universal_counter].legend(loc='best', markerscale=0.25, handlelength=0.25)
                         RS_label = 1

                    else:
                        ax[universal_counter].add_patch(
                            Rectangle((value[0], 0), value[1] - value[0], length[iter][0][2][0][i]['y'], lw=1,
                                      fill=None, color='red'))

            x = vec_values[length[iter][0][1]]['x'][len(vec_values[length[iter][0][1]]['x']) - 1]


            ax[universal_counter].axvline(x, color='cyan', linestyle='--')
            ax[universal_counter].set_ylim([-10, 10])
            ax[universal_counter].set_xlim([0, 30])
            universal_counter += 1
            track_which_ones.append(length[iter][0][1])
            RS_Budget -= 1
         if universal_counter == 19:
             ax[universal_counter].set_xlabel('Computation Time (hr)', fontsize=10, labelpad=5)
     print("tracked stuff",len(track_which_ones))

     plt.savefig("/Users/awahab/Documents/QEERI/IMMI_REVISION/IMMI_NEW_EXPERIMENTS_FIGURES/"+str(k)+".png")

####################################################################################################################### PLOTTING BINS PLOT ############################################################################################################################################
 plt.figure()


 for iter in range(len(length)):
     mean=[]
     std=[]
     fx=[]
     xx=[]
     best = []

     if length[iter][0][1] == plot_reference:
         total_bins=batch_bin(length[iter][0][2][0][len(length[iter][0][2][0]) - 1]['index'], length[iter][0][3])
         print("TOPK-BO #3 data.csv",iter,length[iter][0][2][0])
         print("TOPK-BO #3 data.csv",iter,length[iter][0][2][0][len(length[iter][0][2][0])-1]['index'])
         print("TOPK-BO #3 total bins",total_bins)
         iii=0
         for j in range(length[iter][0][2][0][len(length[iter][0][2][0])-1]['index']):
             y=length[iter][0][2][0][j]['y']
             fx.append(y)
             mean.append(np.mean(fx))
             best.append(min(fx))
             std.append(np.std(fx))
             x=batch_bin(length[iter][0][2][0][j]['index'], length[iter][0][3])
             xx.append(x)
             print("y,x",y,x)
             if iii==0:
                plt.scatter(x,y,color = 'blue',s=10,label="ΔHmix value/s at given batch cycle")
             else:
                 plt.scatter(x, y, color='blue', s=10)
             iii=1

         mean = np.asarray(mean)
         std = np.asarray(std)
         plt.plot(
             xx,
             mean,
             color="grey",label='mean ΔHmix with std dev'

         )
         plt.fill_between(xx, mean + std, mean - std, color="grey", alpha=0.3)

         y_best_value_array=[length[iter][0][2][3][k]['y'] for k in range(len(length[iter][0][2][3]))]
         x_best_value_array=[batch_bin(length[iter][0][2][3][k]['index'],length[iter][0][3]) for k in range(len(length[iter][0][2][3]))]
         print("x_best_value_array",x_best_value_array)
         print("y_best_value_array",y_best_value_array)
         plt.ylabel('ΔHmix (meV/ion)', fontsize=20)
         plt.xlabel('Batch Cycle', fontsize=20)
         plt.plot(xx,best,color='orange',label="Best Value Found")
         if int(str(x_best_value_array[len(x_best_value_array)-1])[-1]) == 1:
             th = "{st}"
         elif int(str(x_best_value_array[len(x_best_value_array)-1])[-1])  == 2:
             th = "{nd}"
         elif int(str(x_best_value_array[len(x_best_value_array)-1])[-1])  == 3:
             th = "{rd}"
         else:
             th = "{th}"
         plt.annotate(
             "Best Value Found {} meV at {}$^{}$ Batch Cycle".format(str(np.round(y_best_value_array[len(y_best_value_array)-1],2)), str(x_best_value_array[len(x_best_value_array)-1]),th),
             xy=(x_best_value_array[len(x_best_value_array)-1] + 0.5, y_best_value_array[len(y_best_value_array)-1]),
             xytext=(x_best_value_array[len(x_best_value_array)-1] + 15, y_best_value_array[len(y_best_value_array)-1]),
             arrowprops=dict(color="green"),
             color="darkgreen",
             bbox=dict(facecolor="white", alpha=1.0),
         )
         plt.legend()


def sigmoid(x):
    return 1 / (1 + np.exp(-x))
def standard_dev(set,mean):
    sum = 0
    for ij in range(len(set)):
        sum += (ij - mean) ** 2
    var = np.sqrt(sum/ len(set))
    return var

def read_data(file,sample_size):
    yvalue=[]
    cut_off=[]
    cutt_off = -1
    min_so_far=[]
    with open(file+'/data.csv', 'r') as f:
     lines= f.readlines()
     print("data lines length", len(lines))
     for i in range(len(lines)):
      if int(lines[i].strip().split(",")[1].split(":")[1]) <= sample_size:
          yvalue.append({'y':float(lines[i].strip().split(",")[0].split(":")[1]),'index':int(lines[i].strip().split(",")[1].split(":")[1])})
          if float(lines[i].strip().split(",")[0].split(":")[1]) <= cutt_off:
              cut_off.append({'y':float(lines[i].strip().split(",")[0].split(":")[1]),'index':int(lines[i].strip().split(",")[1].split(":")[1])})
          if len(min_so_far)==0:
              min_so_far.append({'y':float(lines[i].strip().split(",")[0].split(":")[1]),'index':int(lines[i].strip().split(",")[1].split(":")[1])})
          else:
              if float(lines[i].strip().split(",")[0].split(":")[1]) <= float(min_so_far[(len(min_so_far)-1)]['y']):
                  min_so_far.append({'y':float(lines[i].strip().split(",")[0].split(":")[1]),'index':int(lines[i].strip().split(",")[1].split(":")[1])})
     obtain_min= dict(map(lambda yvalue: (yvalue['y'], yvalue), yvalue))

     print("data lines length adjusted",file,len(yvalue))
     obtain_min = obtain_min[min(obtain_min.keys())]


    return yvalue,obtain_min,cut_off,min_so_far

################ DATA ########################################################################################################
file1="/Users/awahab/Downloads/NEW_HPC_CAL_DATA_2/FDIBO-1/"
file2="/Users/awahab/Downloads/NEW_HPC_CAL_DATA_2/FDIBO-2/"
file3="/Users/awahab/Downloads/NEW_HPC_CAL_DATA_2/FDIBO-3/"
file4="/Users/awahab/Downloads/NEW_HPC_CAL_DATA_2/FDIBO-4/"
file5="/Users/awahab/Downloads/NEW_HPC_CAL_DATA_2/FDIBO-5/"
file6="/Users/awahab/Downloads/NEW_HPC_CAL_DATA_2/TOPK-1/"
file7="/Users/awahab/Downloads/NEW_HPC_CAL_DATA_2/TOPK-2/"
file8="/Users/awahab/Downloads/NEW_HPC_CAL_DATA_2/TOPK-3/"
file9="/Users/awahab/Downloads/NEW_HPC_CAL_DATA_2/TOPK-4/"
file10="/Users/awahab/Downloads/NEW_HPC_CAL_DATA_2/TOPK-5/"
file11="/Users/awahab/Downloads/NEW_HPC_CAL_DATA_2/SBO-1/"
file12="/Users/awahab/Downloads/NEW_HPC_CAL_DATA_2/SBO-2/"


file16="/Users/awahab/Downloads/NEW_HPC_CAL_DATA_2/random_1/"
file17="/Users/awahab/Downloads/NEW_HPC_CAL_DATA_2/random_2/"
file18="/Users/awahab/Downloads/NEW_HPC_CAL_DATA_2/random_3/"
file19="/Users/awahab/Downloads/NEW_HPC_CAL_DATA_2/random_4/"
file20="/Users/awahab/Downloads/NEW_HPC_CAL_DATA_2/random_5/"

file31="/Users/awahab/Downloads/NEW_HPC_CAL_DATA_2/SBO-200SAMPLES-1/"
file32="/Users/awahab/Downloads/NEW_HPC_CAL_DATA_2/SBO-200SAMPLES-2/"
file33="/Users/awahab/Downloads/NEW_HPC_CAL_DATA_2/SBO-200SAMPLES-3/"
file34="/Users/awahab/Downloads/NEW_HPC_CAL_DATA_2/SBO-200SAMPLES-4/"
file35="/Users/awahab/Downloads/NEW_HPC_CAL_DATA_2/SBO-200SAMPLES-5/"
file65="/Users/awahab/Downloads/NEW_HPC_CAL_DATA_2/new/SBO-REF-GBR-LCB1/"
file66="/Users/awahab/Downloads/NEW_HPC_CAL_DATA_2/new/SBO-REF-GBR-LCB4/"
file67="/Users/awahab/Downloads/NEW_HPC_CAL_DATA_2/new/SBO-REF-GPR-LCB/"

sample_size=100

#
p1=[extract_pulses(file1,sample_size),'FDI-BO #1',read_data(file1,sample_size),4]
p2=[extract_pulses(file2,sample_size),'FDI-BO #2',read_data(file2,sample_size),4]
p3=[extract_pulses(file3,sample_size),'FDI-BO #3',read_data(file3,sample_size),4]
p4=[extract_pulses(file4,sample_size),'FDI-BO #4',read_data(file4,sample_size),4]
p5=[extract_pulses(file5,sample_size),'FDI-BO #5',read_data(file5,sample_size),4]
p6=[extract_pulses(file6,sample_size),'TOPK-BO #1',read_data(file6,sample_size),4]
p7=[extract_pulses(file7,sample_size),'TOPK-BO #2',read_data(file7,sample_size),4]
p8=[extract_pulses(file8,sample_size),'TOPK-BO #3',read_data(file8,sample_size),4]
p9=[extract_pulses(file9,sample_size),'TOPK-BO #4',read_data(file9,sample_size),4]
p10=[extract_pulses(file10,sample_size),'TOPK-BO #5',read_data(file10,sample_size),4]
p11=[extract_pulses(file11,sample_size),'S-BO #1',read_data(file11,sample_size),1]
p12=[extract_pulses(file12,sample_size),'S-BO #2',read_data(file12,sample_size),1]
p16=[extract_pulses(file16,sample_size),'RS #1',read_data(file16,sample_size),1]
p17=[extract_pulses(file17,sample_size),'RS #2',read_data(file17,sample_size),1]
p18=[extract_pulses(file18,sample_size),'RS #3',read_data(file18,sample_size),1]
p19=[extract_pulses(file19,sample_size),'RS #4',read_data(file19,sample_size),1]
p20=[extract_pulses(file20,sample_size),'RS #5',read_data(file20,sample_size),1]
p66=[extract_pulses(file66,sample_size),'S-BO #21',read_data(file66,sample_size),1]
p67=[extract_pulses(file67,sample_size),'S-BO #22',read_data(file67,sample_size),1]



###########################################################   IMPROVING STASTICS ##################################################################

plotting_profiler([p1],[p2],[p3],[p4],[p5],[p6],[p7],[p8],[p9],[p10],[p11],[p12],[p16],[p17],[p18],[p19],[p20],[p66],[p67])

###################################################################################################################################################
# qr=[]
# for key in data.keys():
#     print("key",key)
#     qr.append(data[key])
# print(qr)
# plotting_profiler(qr)
plt.show()



