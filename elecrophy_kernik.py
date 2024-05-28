import myokit
from PIL import Image, ImageOps
import os
import glob
import numpy as np
import matplotlib.pyplot as plt

def make_protocol(stim=1,
                  start_time=50,
                  duration=25,
                  period=50,
                  multiplier=3):
    '''Creates protocol of stimulatioms
        Params:
            stim(float) - stimulation amplitude
            start_time(int) - time of first stimulation
            duration(int) - time of each stimulus
            period(int) - period of stimulations
            multiplier(int) - number of periodic stimulations
    '''

    protocol = myokit.Protocol()
    protocol.schedule(stim, start_time, duration, period, multiplier)

    return protocol


def make_simulation(model,
                    protocol, 
                    gx,
                    gy,
                    npixels=(100, 100),
                    step_size=0.01,
                    rl=True
                    ):
    '''
    Insatance simulation class
        Params:
            model(myokit.Model) - model for simulation
            protocol(myokit.Protocol) - stimulation protocol for simulation
            gx(numpy.ndarray) - conductivity along x axis map for simulation
            gy(numpy.ndarray) - conductivity along y axis map for simulation
            npixels(int ot tuple) - size of simulation area
            paced_x(int) - number of x pixels for stimulation
            paced_y(int) - number of x pixels for simulation
            start_x(int) - start x pixel 
            start_y(int) - start y pixel
            step_size(float) - integration time step ms
    '''

    simulation = myokit.SimulationOpenCL(model, protocol, 
                                         ncells=npixels, 
                                         diffusion=True, rl=True,
                                         precision=myokit.DOUBLE_PRECISION,
                                         native_maths=True)
    simulation.set_conductance_field(gx, gy)
    simulation.set_step_size(step_size)

    return simulation

def run_simulation(simulation,
                   protocol,
                   time=50,
                   log_interval=5,
                   log_vars=myokit.LOG_STATE+myokit.LOG_BOUND,
                #    log_vars = ['environment.time','Voltage.Vm'],
                   start_x=0,
                   start_y=0,
                   paced_x = 25,
                   paced_y=100,
                   s1s2=False,
                   s2_paced_x=40,
                   s2_paced_y=200, 
                   after_time=0,
                   circle=False,
                   circle_O=(15,50),
                   circle_R=10):
    '''
    Run simulation 
        Params:
            simulation(myokit.SimulationOpenCL) - simulation to run
            time(int) - time of simulation
            log_interval(int) - step of frames
            log_vars(list of str) - variables to log(refer to mmt file to see possible variables)
            s1s2(bool) - s1s2 protocol simulation or default
            paced_x(int) - number of x pixels for stimulation
            paced_y(int) - number of y pixels for stimulation
            start_x(int) - start x pixel 
            start_y(int) - start y pixel
            s2_paced_x(int) - number of x pixels for s2 stimulation
            s2_paced_y(int) - number of y pixels for s2 stimulation
    '''
    if circle:
        circle_list = [circle_O]
        for i in range(-circle_R, circle_R):
            for j in range(-circle_R, circle_R):
                if (i**2 + j**2) <= circle_R**2:
                    circle_list.append((i+circle_O[0],j+circle_O[1]))
        
        simulation.set_paced_cell_list(circle_list)
    else:
        simulation.set_paced_cells(nx=paced_x, ny=paced_y, x=start_x, y=start_y)
    
    if not s1s2:
        simulation_result = simulation.run(time, log_interval=log_interval, log=log_vars)
        return [simulation_result]
    else:
        period = protocol.head().period()
        t = 0
        simulation_results = []
        
        while t < time:
            simulation_results.append(simulation.run(period, log_interval=log_interval, log=log_vars))
            t += period
            if (t // period) % 2 == 1:
                simulation.set_paced_cells(nx=s2_paced_x, ny=s2_paced_y)
            else:
                simulation.set_paced_cells(nx=paced_x, ny=paced_y, x=start_x, y=start_y)
            period += after_time
        return simulation_results
            

def visual(simulation_result, variable, path='results'):
    if not os.path.isdir(path):
        os.makedirs(path)
    i = 0
    for result in simulation_result:    
        block = result.block2d()
        imgs = block.get2d('Voltage.Vm')
        for array in imgs:
            plt.imsave(f'{path}/{i}.png', array, cmap='gray', vmin=-100, vmax=50)
            i += 1
            
def cond_velocity_2d(simulation_result, pix_size=0.0025, where=999, num_points=10, size=1000, treshold=-30):
    time_range = np.array(simulation_result[0]['environment.time'])
    v_array = np.array([np.array(simulation_result[0][f'{i}.{where}.membrane.V']) for i in range(0, size-1, size // num_points)])
    indecies = np.argmax(v_array >= treshold, axis=1)
    times = time_range[indecies]
    velocities = np.ones(times.shape[0]) * where * pix_size / times * 100
    if len(velocities[velocities > 0]) == 0:
        return 0
    return np.mean(velocities[velocities > 0])

def cond_velocity_2d_hor(simulation_result, pix_size=0.0025, where=999, num_points=10, size=1000, treshold=-30):
    time_range = np.array(simulation_result[0]['environment.time'])
    v_array = np.array([np.array(simulation_result[0][f'{where}.{i}.membrane.V']) for i in range(0, size-1, size // num_points)])
    indecies = np.argmax(v_array >= treshold, axis=1)
    times = time_range[indecies]
    velocities = np.ones(times.shape[0]) * where * pix_size / times * 100
    if len(velocities[velocities > 0]) == 0:
        return 0
    return np.mean(velocities[velocities > 0])
    
def make_gif(path):
    # Create the frames
    frames = []
    imgs = glob.glob(path)
    for i in imgs:
        new_frame = Image.open(os.path.join(path, i))
        frames.append(new_frame)

    # Save into a GIF file that loops forever
    frames[0].save('png_to_gif.gif', format='GIF',
                   append_images=frames[1:],
                   save_all=True,
                   duration=1000, loop=0)
