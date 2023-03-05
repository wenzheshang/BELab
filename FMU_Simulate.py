from pyfmi import load_fmu

import numpy as np
import matplotlib.pyplot as plt


#可以把这里做成一个class的形式，进行调用，或者干脆也集成到main里面


Tstart = 0 # The start time.
Tend = 900

model = load_fmu('FMU/Buildings_Vadilation_vadilation.fmu')
result_dict = model.get_model_variables()
result = list(result_dict.keys())
print(result[0:10])
ms =  open(('FMU/variable_list.txt'),'w')

model.set(['TFloor.T'], [303.15])
model.setup_experiment(start_time = Tstart) # Set the start time to Tstart
model.enter_initialization_mode()
model.exit_initialization_mode()
eInfo = model.get_event_info()
eInfo.newDiscreteStatesNeeded = True
#Event iteration
# while eInfo.newDiscreteStatesNeeded == True:
#     model.enter_event_mode()
#     model.event_update()
#     eInfo = model.get_event_info() #这一部分似乎没啥用处

model.enter_continuous_time_mode()

# Get Continuous States
x = model.continuous_states
# Get the Nominal Values
x_nominal = model.nominal_continuous_states
# Get the Event Indicators
event_ind = model.get_event_indicators()

# Values for the solution
# Retrieve the valureferences for the values 'CFD_roo.Room_MeanT'

vref = [model.get_variable_valueref('CFD_roo.Room_MeanT')]+ \
     [model.get_variable_valueref('SupplyAir.T_in')]
t_sol = [Tstart]
sol = [model.get_real(vref)]
a = model.get_real(vref)#[0]
time = Tstart
Tnext = Tend # Used for time events
dt = 0.01
#value = 283.15

number = 0

while time < Tend and not model.get_event_info().terminateSimulation:
    #Compute the derivative of the previous step f(x(n), t(n))
    dx = model.get_derivatives()
    
    # Advance
    h = min(dt, Tnext-time)
    time = time + h
    
    # Set the time
    model.time = time
    
    # Set the inputs at the current time (if any)
    #model.set(vref, value)
    
    # Set the states at t = time (Perform the step using x(n+1)=x(n)+hf(x(n), t(n))
    x = x + h*dx 
    model.continuous_states = x

    # Get the event indicators at t = time
    event_ind_new = model.get_event_indicators()
    
    # Inform the model about an accepted step and check for step events
    step_event = model.completed_integrator_step()
    
    # Check for time and state events
    time_event = abs(time-Tnext) <= 1.e-10
    state_event = True if True in ((event_ind_new>0.0) != (event_ind>0.0)) else False

    # Event handling
    if step_event or time_event or state_event:
        model.enter_event_mode()
        eInfo = model.get_event_info()
        eInfo.newDiscreteStatesNeeded = True
        model.set(set_variable, set_value)
        if number >= communicate_time:
            iii = -iii
            number = 0
        number = number + 0.01
        
        # Event iteration
        while eInfo.newDiscreteStatesNeeded:
            model.event_update('0') # Stops at each event iteration
            eInfo = model.get_event_info()
        # Retrieve solutions (if needed)
    if eInfo.newDiscreteStatesNeeded:
        pass
    
    # Check if the event affected the state values and if so sets them
    if eInfo.valuesOfContinuousStatesChanged:
        x = model.continuous_states
    
    # Get new nominal values.
    # if eInfo.nominalsOfContinuousStatesChanged:
    #     atol = 0.01*rtol*model.nominal_continuous_states
    
    # Check for new time event
    if eInfo.nextEventTimeDefined:
        Tnext = min(eInfo.nextEventTime, Tend)
    else:
        Tnext = Tend
    model.enter_continuous_time_mode()

    event_ind = event_ind_new

    # Retrieve solutions at t=time for outputs
    # bouncing_fmu.get_real,get_integer,get_boolean,get_string (valueref)

    t_sol += [time]
    sol += [model.get_real(vref)]


# model.get_model_variables()
# res = model.simulate(final_time=720)
# t = res['time']
# x1 = res['SupplyAir.T_in']
plt.subplot(211)
plt.plot(t_sol,np.array(sol)[:,0])
plt.subplot(212)
plt.plot(t_sol,np.array(sol)[:,1])
plt.savefig('/FMU/FMU_result.svg')


# dir_result = 'F:\\Thinking\\modelica+CFD\\Vadilation'
# result_path = os.path.join(dir_result,'vadilation.mat')
# r = Reader(result_path,'dymola')
# (time,RoomMeanT) = r.values('CFD_roo.air.heaPorAir.T')
# (time_2,SupplyT) = r.values('SupplyAir.T_in')
# mydataframe = pd.DataFrame({'time_roo':time, 'RoomMeanT':RoomMeanT, 'time_Sup':time_2, 'SupplyT': SupplyT})
# mydataframe.to_csv(os.path.join(dir_result, 'data/result.csv'))
# signal_simu_time = [1]
# dymola = DymolaInterface()
# smash_signal = [0]


# def simulate(value, name):

#     dir_result = 'F:\\Thinking\\modelica+CFD\\Vadilation\\data'
#     result_path = os.path.join(dir_result,'vadilation.mat')
#     modelicaPath = pathlib.Path(os.environ["DymolaPath"])
#     #Library import
#     dirBuilding = os.path.join(modelicaPath,"Modelica/Library/Buildings-v8.0.0/Buildings 8.0.0")
#     #open Library
#     dymola.openModel(path=os.path.join(dirBuilding, 'package.mo'))
#     dymola.openModel(path='F:\\Thinking\\modelica+CFD\\Vadilation\\Vadilation.mo')
#     problemName = 'Vadilation.vadilation'#'Plant.plant_rectify_0105_correct'
#     ResultValue = []
#     demo_name = 'demo_results'+str(signal_simu_time[0])
#     (dymola_setName, dymola_setValue) = fluent_to_dymola(value, name)
#     print(dymola_setName, dymola_setValue)
#     try:
#         result = dymola.simulateExtendedModel(
#                     problem= problemName,
#                     startTime=0,
#                     stopTime=180,
#                     numberOfIntervals=0,
#                     outputInterval=0.0,
#                     method="Dassl",
#                     tolerance=0.0001,
#                     fixedstepsize=0.0,
#                     resultFile=os.path.join(dir_result,demo_name),
#                     initialNames=dymola_setName,
#                     initialValues=dymola_setValue,
#                     autoLoad=True
#                     )
#     except:
#         print('error')
#         log = dymola.getLastError()
#         f =  open(os.path.join(dir_result,'error.txt'),'w')
#         f.write(log)
#         f.close()
#         return
#     result_path = os.path.join(dir_result,demo_name+'.mat')
#     r = Reader(result_path,'dymola')

#     result_name = 'reslut'+str(signal_simu_time[0])
#     signal_simu_time[0] = signal_simu_time[0] + 1

#     ResultVarName = r.varNames() #获取所有结果变量名
#     for i in range(len(ResultVarName)):
#         (t,r_ser) = r.values(ResultVarName[i])
#         ResultValue.append(r_ser[-1])
#     mydataframe = pd.DataFrame({'VarName':ResultVarName,'Value':ResultValue})
#     mydataframe.to_csv(os.path.join(dir_result, result_name+'.csv'))#将所有结果保存到.csv文件中，以备下次读取

# def fluent_to_dymola(value, name):
#     Exchange_data_value = [value]
#     Exchange_data_name = [name]
#     return Exchange_data_name, Exchange_data_value

# for i in range(5):
#     value = 283.15+10*i
#     name = 'CFD_roo.Room_MeanT'
#     simulate(value, name)



# dymola = DymolaInterface()
# ResultValue = []

# result_path = 'C:\\Users\\Administrator\\Desktop\\result.csv'
# Exchange_data_name = []
# Exchange_data_value = []
# with open(result_path) as csvfile:
#     reader = csv.DictReader(csvfile)
#     for row in reader:      
#         Exchange_data_name.append(row['VarName'])
#         Exchange_data_value.append(float(row['Value']))#将csv中结果读入列表

# modelicaPath = pathlib.Path(os.environ["DymolaPath"])
# #Library import
# dirBuilding = os.path.join(modelicaPath,"Modelica/Library/Buildings-v8.0.0/Buildings 8.0.0")
# #open Library
# dymola.openModel(path=os.path.join(dirBuilding, 'package.mo'))
# dymola.openModel(path='F:/Thinking/modelica+CFD/Plant.mo')
# problemName = 'Plant.plant_rectify_0105_correct'#'vadilation'


# dymola_setName = Exchange_data_name
# dymola_setValue = Exchange_data_value
# endT = 1800
# intervals = 10
# step_time = endT/intervals


# result = dymola.simulateExtendedModel(
#     problem= problemName,
#     startTime=0,
#     stopTime=step_time,
#     numberOfIntervals=0,
#     outputInterval=0.0,
#     method="Dassl",
#     tolerance=0.0001,
#     fixedstepsize=0.0,
#     resultFile=os.path.join('/test_result','demo_results'),
#     initialNames=dymola_setName,
#     initialValues=dymola_setValue,
#     autoLoad=True
#     )

# status = result[0]
# if not status:
#     print('error')
#     log = dymola.getLastError()
#     f =  open(os.path.join('C:\\Users\\Administrator\\Desktop\\error.txt'),'w')
#     f.write(log)
#     f.close()

# else:
#     #成功模拟后输出结果部分,加保存excel功能
#     #以下代码保存excel文件
#     result_path = os.path.join('/test_result','demo_results.mat')
#     r = Reader(result_path,'dymola')

#     result_name = 'reslut.csv'
    

#     ResultVarName = r.varNames() #获取所有结果变量名
#     for i in range(len(ResultVarName)):
#         (t,r_ser) = r.values(ResultVarName[i])
#         ResultValue.append(r_ser[-1])
#     mydataframe = pd.DataFrame({'VarName':ResultVarName,'Value':ResultValue})
#     mydataframe.to_csv(os.path.join('/test_result', result_name))#将所有结果保存到.csv文件中，以备下次读取