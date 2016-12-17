import math
import time
import random

class Distribution:
    def __init__(self, type='gaussian'):
        print 'initializing %s distribution handler'%(type)
        self.type = type
        if type == 'gaussian':
            self.initial_constants = [476,0.54,223.6]
            self.arg_cnt = 2
        elif type == 'cauchy':
            self.initial_constants = [0.05,6.56,1.97,0.27]
            self.arg_cnt = 4
        elif type == 'double gaussian':
            self.initial_constants = [787.58,5.01,822.81,152.29,-12.49,249.53]
            self.arg_cnt = 4
    def randomize_initial_constants(self, seed='None'):
        if seed=='None':
            seed = time.clock()
            random.seed(seed)
        else:
            random.seed(seed)
        constants = []
        if self.type == 'gaussian':
            constants.append([300.0, 600.0])
            constants.append([0.1,   0.9])
            constants.append([200.0, 300.0])
        elif self.type == 'cauchy':
            constants.append([0.00,2.00])
            constants.append([0.00,10.00])
            constants.append([0.00,5.00])
            constants.append([0.00,2.00])
        elif self.type == 'double gaussian':
            constants.append([600,900])
            constants.append([0.05,0.25])
            constants.append([400.0,600.0])
            constants.append([100,300])
            constants.append([0.05,0.15])
            constants.append([100,200])
        for i in range(0,self.arg_cnt):
            self.initial_constants[i] = random.uniform(constants[i][0],constants[i][1])
            
    def get_parameter_approximations(self, panel, constants='None'):
        # peaks is a list of indices to panels .x and .y for peaks recognized
        if constants=='None':
            args_in = self.initial_constants
        else:
            args_in = constants
        args_out = []
        for peak in panel.peaks:
            x = panel.x[peak]
            y = panel.y[peak]
            tmp_arg = self.get_parameter_approximation(panel.x[peak], panel.y[peak], panel.x[-1], args_in)
            args_out.append(tmp_arg)
        return args_out

    def get_parameter_approximation(self, x, y, last_x, constants='None'):
        if constants == 'None':
            args_in = self.initial_constants
        else:
            args_in = constants
        tmp_arg = []
        if self.type == 'gaussian':
            tmp_arg.append(x/float(args_in[0]))     # the variance of peaks increases with mass
            tmp_arg.append(float(args_in[1])*y + (y * x / float(args_in[2])))       # height increases with peak height and mass
        elif self.type == 'cauchy':
            tmp_arg.append(args_in[0]*y)
            tmp_arg.append(args_in[1])
            tmp_arg.append(args_in[2])
            tmp_arg.append(args_in[3])
        elif self.type == 'double gaussian':
            tmp_arg.append((x/float(args_in[0])) )     # the variance of peaks increases with mass
            tmp_arg.append(float(args_in[1])*y + (y * (x**(1.0-(0.1*x/last_x))) / float(args_in[2])))       # height increases with peak height and mass
            tmp_arg.append((x/float(args_in[3])) )     # the variance of peaks increases with mass 
            tmp_arg.append(float(args_in[4])*y + (y * (x**(1.0+(0.1*x/last_x))) / float(args_in[5])))       # height increases with peak height and mass
        return tmp_arg
    
    def calculate_initial_distributions(self, parms, panel, peak_indices):
        self.calculate_distributions(parms, panel, peak_indices, 0, len(panel.x), 1)
    
    def calculate_distributions(self, parms, panel, peak_indices, start, end, initialize=0):
        # for each coordinate, calculate each peak's contribution, save sum in dist_y
        cnt = 0
        if initialize:   #doing a full one... create dist_y from scratch
            # initialize the distribution sum array
            panel.dist_y = []
            for coor_index in range(start, end):
                ht = 0.0
                cnt = 0
                for peak_index in peak_indices:
                    x = panel.x[peak_index]
                    if abs(panel.x[coor_index] - x) < 60.0*(x/500.0):
                        ht = ht + self.get_distribution_height(parms[cnt], panel, coor_index, peak_index)
                    cnt = cnt + 1
                panel.dist_y.append(ht)
        else:
            # don't initialize
            for coor_index in range(start, end):
                ht = 0.0
                cnt = 0
                for peak_index in peak_indices:
                    x = panel.x[peak_index]
                    if abs(panel.x[coor_index] - x) < 60.0*(x/500.0):
                        ht = ht + self.get_distribution_height(parms[cnt], panel, coor_index, peak_index)
                    cnt = cnt + 1
                panel.dist_y[coor_index] = ht
                
    def get_distribution_height(self, args, panel, coor, peak):
        if self.type == 'gaussian':
            return args[1] * (1.0/(math.sqrt(2*math.pi)*args[0])) * math.exp(-(((panel.x[coor]/args[0])- ((panel.x[peak])/args[0]))**2)/2.0)
        elif self.type == 'cauchy':
            return ((abs(1/args[0])**3)/(math.pi*((abs(args[3]/args[0])**args[1]) + ((abs(panel.x[coor]-panel.x[peak]))**args[2]))))
        elif self.type == 'double gaussian':
            return (args[1] * (1.0/(math.sqrt(2*math.pi)*args[0])) * math.exp(-(((panel.x[coor]-panel.x[peak])/args[0])**2)/2.0)) + (args[3] * (1.0/(math.sqrt(2*math.pi)*args[2])) * math.exp(-(((panel.x[coor] - panel.x[peak])/args[2])**2)/2.0))
        
    def plot_distributions(self, panel):
        print "plotting %d %s distributions"%(len(panel.peak_args), self.type)
        panel.plot.delete('peak_plot')
        for peak_index in range(0,len(panel.peak_args)):
            print '.',
            for coor_index in range(len(panel.x)-1):
                x = panel.x[panel.peaks[peak_index]]
                if abs(panel.x[coor_index] - x) < 50.0*(x/700.0):
                    ht      = self.get_distribution_height(panel.peak_args[peak_index], panel, coor_index, panel.peaks[peak_index])
                    next_ht = self.get_distribution_height(panel.peak_args[peak_index], panel, coor_index+1, panel.peaks[peak_index])
                    if ht > 0.1 or next_ht > 0.1:
                        panel.plot.create_line(panel.x[coor_index], ht, panel.x[coor_index+1], next_ht, width=1, fill='blue', tags='peak_plot')
        print 
        # scale the image
        scaleval = panel.parent.parent.heightScaler.get()
        panel.plot.scale("peak_plot", 0, 0, 1, -1.5*scaleval)
        panel.plot.move("peak_plot", 0, (panel.height * 0.80))
            
            