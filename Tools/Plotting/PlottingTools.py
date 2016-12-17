# python imports
import parms
import string
import sys
import math
import time
import os.path
import tkFileDialog
import pickle
# internal imports
sys.path.append('./Tools/Math')
import Distribution
# dependency imports
sys.path.append('./Dependencies')   # just in case its not built in
import scipy.stats
import random

# These tools perform operations on sets of .csv files, in the form of class Plot
# originally written for mass spectrometry, the methods are general to peak
# recognition and quantification. The default parameters are set for mass spec though.

def scale_equivalent_observations_by_average_ratio(plots_in, resolution=0.005):
    for p in plots_in:
        p.recognize_peaks(800)
    # locate the peaks that are present in the first spectrum
    weights_present = {}
    p = plots_in[0]
    for ind in p.peaks:
        weight = p.x[ind]
        weights_present[weight] = [0]
        
    for i in range(1,len(plots_in)):
        for ind in plots_in[i].peaks:
            keys = weights_present.keys()
            weight = plots_in[i].x[ind]
            for key in keys:
                if key-weight < weight*resolution:
                    if i not in weights_present[key]:
                        weights_present[key].append(i)
                        
    keys = weights_present.keys()
    keys.sort()
    delkeys = []
    for key in keys:
        for i in range(len(plots_in)):
            if i not in weights_present[key]:
                if key not in delkeys:
                    delkeys.append(key)
                
    for key in delkeys:
        del weights_present[key]
        
    keys = weights_present.keys()
    keys.sort()
    
    # use these keys to get indices for the first plot,
    p = plots_in[0]
    first_indices = []
    for key in keys:
        for i in range(len(p.x)):
            if p.x[i] >= key:
                first_indices.append(i)
                break

    # then use these indices with the get_trough_indices_for_peak_index function to get indices
    indices_table = []
    for i in range(len(first_indices)):
        indices_table.append([])
        for p in plots_in:
            a,b,c = p.get_trough_indices_for_peak_index(first_indices[i])
            indices_table[-1].append(b)

    f = open('./ratios.txt', 'a')
    # scale each observation to the first.
    for i in range(1,len(plots_in)):
        ratio_sum = 0.0
        # use the indices to get the average ratio
        lines = ['ratios plot %s\n'%(i), '']
        for j in range(len(indices_table)):
            ratio = plots_in[0].y[indices_table[j][0]] / plots_in[i].y[indices_table[j][i]]
            lines[-1] += '%5.3f, '%(ratio)
            ratio_sum += ratio
        lines[-1] += '\n'
        ratio = 100.0*(ratio_sum / len(indices_table))
        plots_in[i].hard_rescale_y(ratio)
    f.writelines(lines)
    f.close()

    # now calculate the average peak height for the set of consensus peaks and print
    for i in range(len(plots_in)):
        sum = 0.0
        for j in range(len(indices_table)):
            sum += plots_in[i].y[indices_table[j][i]]

def scale_equivalent_observations_by_height_sum(plots_in, resolution=0.005):
    for p in plots_in:
        p.recognize_peaks(800)
    
    weights_present = {}
    p = plots_in[0]
    for ind in p.peaks:
        weight = p.x[ind]
        weights_present[weight] = [0]
        
    for i in range(1,len(plots_in)):
        for ind in plots_in[i].peaks:
            keys = weights_present.keys()
            weight = plots_in[i].x[ind]
            for key in keys:
                if key-weight < weight*resolution:
                    if i not in weights_present[key]:
                        weights_present[key].append(i)
                        
    keys = weights_present.keys()
    keys.sort()
    delkeys = []
    for key in keys:
        for i in range(len(plots_in)):
            if i not in weights_present[key]:
                if key not in delkeys:
                    delkeys.append(key)
                
    for key in delkeys:
        del weights_present[key]
        
    keys = weights_present.keys()
    keys.sort()
    
    # use these keys to get indices for the first plot,
    p = plots_in[0]
    first_indices = []
    for key in keys:
        for i in range(len(p.x)):
            if p.x[i] >= key:
                first_indices.append(i)
                break

    # then use these indices with the get_trough_indices_for_peak_index function to get indices
    # the following peaks are being used in the initial scaling
    indices_table = []
    row_weights = []
    for i in range(len(first_indices)):
        indices_table.append([])
        row_weights.append(plots_in[0].x[first_indices[i]])
        for p in plots_in:
            a,b,c = p.get_trough_indices_for_peak_index(first_indices[i])
            indices_table[-1].append(b)
            
    if len(indices_table) == 0:
        print 'no consensus set of peaks to use in height sum. aborting initial scaling'
        return

    if len(indices_table) > 1:
        row_avgs = []
        row_sds = []
        for i in range(len(indices_table)):
            row_avgs.append(0.0)
            row_sds.append(0.0)

        for i in range(len(indices_table)):
            for j in range(len(plots_in)):
                row_avgs[i] += plots_in[j].y[indices_table[i][j]]/float(len(indices_table[i]))

        for i in range(len(indices_table)):
            for j in range(len(plots_in)):
                row_sds[i] += (plots_in[j].y[indices_table[i][j]]-row_avgs[i])**2
            
        for i in range(len(row_sds)):
            row_sds[i] = math.sqrt(row_sds[i]/float(len(row_sds)-1))

        print 'pre-normalization analysis for scaling by height sum'
        for i in range(len(row_weights)):
            print '%6.1f - '%(row_weights[i]),
            for j in range(len(indices_table[i])):
                print '%5.1f '%(plots_in[j].y[indices_table[i][j]]),
            print ' avg=%5.1f sd=%5.1f val=%5.1f'%(row_avgs[i], row_sds[i], row_sds[i]/row_avgs[i])

        sums = []
        lsms = len(indices_table[0])
        for j in range(len(indices_table[0])):
            sums.append(0.0)
        for j in range(lsms):
            for item in indices_table:
                sums[j] += plots_in[j].y[item[j]]
            
        column_sums = []
        for sum in sums:
            column_sums.append(sum)
        print "         ",
        for avg in column_sums:
            print '%5.1f '%(avg),
        print

    # done pre-normalization analysis and printout
    
    # collect the height sums
    sums = []
    for p in plots_in:
        inisum = 0.0
        for j in range(len(indices_table)):
            inisum += p.y[indices_table[j][0]]
        sums.append(inisum)

    # find the largest
    largest = None
    size = 0
    for i in range(len(sums)):
        if sums[i] > size:
            size = sums[i]
            largest = i
    
    # scale each observation to the first.
    inisum = 0.0
    for j in range(len(indices_table)):
        inisum += plots_in[largest].y[indices_table[j][largest]]
    
    for i in range(0,len(plots_in)):
        if i != largest:
            second_sum = 0.0
            for j in range(len(indices_table)):
                second_sum += plots_in[i].y[indices_table[j][i]]
            if second_sum > 0.0:
                ratio = 100.0*( inisum / second_sum )
                plots_in[i].hard_rescale_y(ratio)
            else:
                print 'no peaks found, no rescale to apply'

    # now calculate the average peak height for the set of consensus peaks and print
    for i in range(len(plots_in)):
        sum = 0.0
        for j in range(len(indices_table)):
            sum += plots_in[i].y[indices_table[j][i]]
        print 'plot %s, peakheight sum = %s'%(plots_in[i].label, sum)

    # begin post-normalization analysis and printout


    if len(indices_table) > 1:
        row_avgs = []
        row_sds = []
        for i in range(len(indices_table)):
            row_avgs.append(0.0)
            row_sds.append(0.0)

        for i in range(len(indices_table)):
            for j in range(len(plots_in)):
                row_avgs[i] += plots_in[j].y[indices_table[i][j]]/float(len(indices_table[i]))

        for i in range(len(indices_table)):
            for j in range(len(plots_in)):
                row_sds[i] += (plots_in[j].y[indices_table[i][j]]-row_avgs[i])**2
            
        for i in range(len(row_sds)):
            row_sds[i] = math.sqrt(row_sds[i]/float(len(row_sds)-1))

        print 'post-normalization analysis'
        for i in range(len(row_weights)):
            print '%6.1f - '%(row_weights[i]),
            for j in range(len(indices_table[i])):
                print '%5.1f '%(plots_in[j].y[indices_table[i][j]]),
            print ' avg=%5.1f sd=%5.1f val=%5.1f'%(row_avgs[i], row_sds[i], row_sds[i]/row_avgs[i])

        sums = []
        lsms = len(indices_table[0])
        for j in range(len(indices_table[0])):
            sums.append(0.0)
        for j in range(lsms):
            for item in indices_table:
                sums[j] += plots_in[j].y[item[j]]
        column_sums = []
        for sum in sums:
            column_sums.append(sum)
        print "         ",
        for avg in column_sums:
            print '%5.1f '%(avg),
        print
    

def normalize_plots_by_unreactive_fragments(plots_in, possible_fragments, resolution = 0.002):
    # recognize peaks in the spectra
    for p in plots_in:
        p.recognize_peaks(800)

    # figure out which possible_fragments correspond to the recognized peaks
    associations = []
    i = 0
    for p in plots_in:
        associations.append([])
        for peak in p.peaks:
            weight = p.x[peak]
            best_frag = None
            closest_distance=100000
            for pf in possible_fragments:   # a set of weights
                d = abs(weight-pf)
                if d < resolution*pf:
                    if d < closest_distance:
                        closest_distance = d
                        best_frag = pf
            if best_frag:
                associations[-1].append([weight, best_frag, peak])
        i += 1

    # figure out which ones are present in any of the spectra
    all_weights = {}
    for set in associations:
        redundancy_monitor = []
        for a in set:
            frag = a[1]
            if frag not in redundancy_monitor:
                redundancy_monitor.append(frag)
                if frag not in all_weights.keys():
                    all_weights[frag] = 1
                else:
                    all_weights[frag] += 1
    

    keys = all_weights.keys()
    keys.sort()

    # print out the heights of fragments in associations
    for target_weight in keys:
        if all_weights[target_weight] == len(plots_in):     # if the association is present for all spectra
            print '%s present in %s spectra'%(target_weight, all_weights[target_weight])

    heights_table = []

    for target_weight in keys:
        heights_table.append([])
        if all_weights[target_weight] == len(plots_in):     # if the association is present for all spectra
            i = 0
            for set in associations:                # one for each spectrum
                for a in set:                       # for each of its peak/fragment associations
                    if a[1] == target_weight:       # if the association involves this fragment
                        heights_table[-1].append(plots_in[i].y[a[2]])
                        break
                i += 1

    f = open('./rxn_data.txt', 'wb')
    pickle.dump(heights_table, f)
    f.close()

def calibrate_x_spread_by_clustering(plot_in, possible_fragments, resolution = 0.005):
    # the following function searches for peak/fragment associations wihin range factor*resolution using shift
    def get_associations(factor, shift):
        associations = []
        for peak in plot_in.peaks:
            weight = plot_in.x[peak] + shift
            height = plot_in.y[peak]
            if height < 1.0:
                continue
            best_frag = None
            closest_distance=100000
            for pf in possible_fragments:   # a set of weights
                d = abs(weight-pf)
                if d < factor*resolution*pf:
                    if d < closest_distance:
                        closest_distance = d
                        best_frag = pf
            if best_frag:
                associations.append([weight-shift, best_frag])
        return associations
        
    # recognize peaks in the spectra
    plot_in.recognize_peaks(800)
    initial_association_count = len(get_associations(1.0, 0.0))
    associations  = get_associations(1.5, 0.0)
    print '%s seconds -- %2s associations at 1.5, %2s associations at 1.0'%(plot_in.label, len(associations), initial_association_count)
    avg_shift = 0
    if len(associations) == 1:
        avg_shift = associations[0][1]-associations[0][0]
    elif len(associations) > 1:
        clusters = {}
        # try to use each distance as a new possible shift
        for i in range(len(associations)):
            shift = associations[i][1]-associations[i][0]   # frag wt - peak wt
            j = 0
            clusters[shift] = get_associations(0.2, shift)

        keys = clusters.keys()
        keys.sort()
        #for shift in keys:
        #    print 'shifting %4.1f gives %s associations'%(shift, len(clusters[shift]))
    
        # find the largest cluster
        best_key = -1
        best_count = 0
        for key in keys:
            if len(clusters[key]) > best_count:
                best_key = key
                best_count = len(clusters[key])

        # figure out how many clusters have this amount
        ties = 0
        for key in keys:
            if len(clusters[key]) == best_count:
                ties += 1

        # if there's a tie, find the shift that's closest to zero
        if ties > 0:
            smallest_key_index = None
            for i in range(0,len(keys)):
                if len(clusters[keys[i]]) == best_count:   # if its one of the ties, consider it
                    if smallest_key_index:
                        if abs(keys[i]) < abs(keys[smallest_key_index]):    # if closer to zero
                            smallest_key_index = i
                    else:
                        smallest_key_index = i  # initialize
            best_key = keys[smallest_key_index]

        associations  = get_associations(0.2, best_key)
        
        avg_shift = 0.0
        for association in associations:
            avg_shift += association[1] - association[0]
        avg_shift /= len(associations)
                
        new_association_count = len(get_associations(1.0, avg_shift))

        print '%2s 1.0, %2s 0.2 -- shift %s, avg shift %s'%(new_association_count, len(associations), best_key, avg_shift)

    return avg_shift
        

def calibrate_x_spread_to_possible_fragments(plot_in, possible_fragments, resolution = 0.005):
    # recognize peaks in the spectra
    plot_in.recognize_peaks(800)
    # figure out which possible_fragments correspond to the recognized peaks
    associations = []
    for peak in plot_in.peaks:
        weight = plot_in.x[peak]
        height = plot_in.y[peak]
        if height < 3.0:
            continue
        best_frag = None
        closest_distance=100000
        for pf in possible_fragments:   # a set of weights
            d = abs(weight-pf)
            if d < resolution*pf:
                if d < closest_distance:
                    closest_distance = d
                    best_frag = pf
        if best_frag:
            associations.append([weight, best_frag])

    start_distance = 0.0
    for ass in associations:
        start_distance += ass[1]-ass[0]
    if len(associations) != 0:
        start_distance /= len(associations)
    else:
        print 'no peak/fragment associations to work with... aborting calibration'
        return
        
    print '%3s peak/fragment associations from %3s peaks for calibration - '%(len(associations), len(plot_in.peaks)),

    best_so_far = 10000
    best_args   = []
    best_ini    = []

    if len(associations) == 0:
        print 'no peak-fragment associations; do not calibrate'
        return
    elif len(associations) == 1:
        print 'only one peak-fragment association. calibrate'
        b = associations[0][0] - associations[0][1]
    else:
        for i in range(10):
            # get the best-fitting y = mx+b equation
            seed = time.clock()
            random.seed(seed)
            initial_b = random.uniform(start_distance-(0.2*start_distance), start_distance+(0.2*start_distance))
        
            c_lambda = lambda args: _get_b_fitness(args, associations)
            opt_args = scipy.stats.optimize.fmin(c_lambda, [initial_b], args=(), xtol=1e-2, ftol=1e-2, maxiter=1000, full_output=0, disp=0)
            end_sd = c_lambda([opt_args])
            if end_sd < best_so_far:
                best_so_far = end_sd
                best_args = opt_args

        print 'best %s  sd %s'%(best_args, best_so_far)

        # modify the x-values accordingly
        b = best_args[0]

    new_xvals = []
    for x in plot_in.x:
        new_xvals.append(x + b)
    plot_in.x = new_xvals
    
def _get_b_fitness(args, blist):
    b = args[0]
    if len(blist) == 1:
        return blist[0][1]-blist[0][0] + b
    
    distances = []
    total_dist = 0.0
    for xy in blist:
        val = xy[1]-(xy[0]+b)
        distances.append( abs(val) )
        total_dist += val
        
    avg_d = total_dist / float(len(blist))

    # get the SD
    sd = 0.0
    for d in distances:
        sd += (d-avg_d)**2
    ret_val = math.sqrt(sd/(len(distances)-1))
    return ret_val
    
def _get_mb_fitness(args, mblist):
    m = args[0]
    b = args[1]
    distances = []
    total_dist = 0.0
    for xy in mblist:
        x = xy[0]
        y = xy[1]
        val = abs(y-((m*x)+b))
        distances.append( val )
        total_dist += val
        
    avg_d = total_dist / float(len(mblist))

    # get the SD
    sd = 0.0
    for d in distances:
        sd += (d-avg_d)**2
    ret_val = math.sqrt(sd/(len(distances)-1))
    return ret_val

def max_of_plots(plots_in, filename_out=None):
    """ first interpolates all plots, then saves and returns the max from each point """
    newP = Plot()
    for i in range(len(plots_in[0].x)):
        newP.x.append(0.0)
        newP.y.append(0.0)
    x = plots_in[0].x
    numplots = len(plots_in)
    for i in range(len(x)):
        maxy = 0.0
        for p in plots_in:
            if p.y[i] > maxy:
                maxy = p.y[i]
        newP.x[i] = x[i]
        newP.y[i] = maxy
    if filename_out:
        if len(filename_out) == 0:
            filename_out = 'average_spectra'
        tokens = string.split(filename_out, '.')
        if tokens[-1] != 'csv':
            filename_out = filename_out + '.csv'
        newP.source_file = filename_out
        newP.write(filename_out)
        return newP
    else:
        return newP

def average_plots(plots_in, filename_out=None):
    """ first interpolates all three, then calculates and saves an average Plot """
    orig_length = len(plots_in[0].x)
    for plot in plots_in:
        if len(plot.x) != orig_length:
            print 'interpolating'
            interpolate_plots(plots_in)
        break
        
    newP = Plot()
    for i in range(len(plots_in[0].x)):
        newP.x.append(0.0)
        newP.y.append(0.0)
    x = plots_in[0].x
    numplots = len(plots_in)
    for i in range(len(x)):
        avey = 0.0
        for p in plots_in:
            avey += p.y[i]
        newP.x[i] = x[i]
        newP.y[i] = avey/numplots
    if filename_out:
        if len(filename_out) == 0:
            filename_out = 'average_spectra'
        tokens = string.split(filename_out, '.')
        if tokens[-1] != 'csv':
            filename_out = filename_out + '.csv'
        newP.source_file = filename_out
        newP.write(filename_out)
        return newP
    else:
        return newP
                    
def interpolate_plots(plots_in):
    # first interpolate the later panels to the first
    # use a simple triangles method for now, but the vtk code can interpolate Bezier curves
    # first find the latest start and earliest end
    latest_start = -10000.0
    earliest_end = 1000000.0
    for p in plots_in:
        if p.x[0] > latest_start:
            latest_start = p.x[0]
        if p.x[-1] < earliest_end:
            earliest_end = p.x[-1]
    # now delimit the first mspanel to these bounds
    p1 = plots_in[0]
    start_ind = 0
    end_ind = len(p1.x)
    for x_ind in range(len(p1.x)):
        if p1.x[x_ind] >= latest_start:
            start_ind = x_ind
            break
    for x_ind in range(len(p1.x)-1, -1, -1):
        if p1.x[x_ind] <= earliest_end:
            end_ind = x_ind
            break
    p1.x = p1.x[start_ind:end_ind+1]
    p1.y = p1.y[start_ind:end_ind+1]
    for p2 in plots_in[1:]:
        new_xs = []
        new_ys = []
        last_2ind = 1
        for x1 in p1.x:
            # find the nearest pair of coordinates in p2
            for x2_ind in range(last_2ind, len(p2.x)):
                if p2.x[x2_ind] >= x1:
                    # use triangles to interpolate
                    ratio = (x1-p2.x[x2_ind-1]) / (p2.x[x2_ind]-p2.x[x2_ind-1])
                    new_xs.append(x1)
                    new_ys.append(p2.y[x2_ind-1] + ratio * (p2.y[x2_ind]-p2.y[x2_ind-1]))
                    # this should work
                    if p2.x[x2_ind] > x1:
                        last_2ind = x2_ind-1
                    else:
                        last_2ind = x2_ind
                    break
        p2.x = new_xs
        p2.y = new_ys
        print '.',
        
def subtract_plot(plots_in, plot_to_subtract):
    for plot in plots_in:
        if len(plot.y) != len(plot_to_subtract.y):
            plots_in.append(plot)
            interpolate_plots(plots_in)
            plot_to_subtract = plots_in[-1]
            plots_in = plots_in[:-1]
    for plot in plots_in:
        for i in range(len(plot_to_subtract.y)):
            plot.y[i] = plot.y[i] - plot_to_subtract.y[i]
            
def normalize_and_subtract_plot(plots_in, plot_to_subtract):
    """ to subtract the controls,
            for each control,
            1. recognize peaks in the control
            2. find the biggest height ratio between the plots at any of these peaks
            3. normalize the plot_to_subtract by this ratio
            4. subtract
    """
    plot_to_subtract.recognize_peaks(1000)
    if len(plot_to_subtract.peaks) == 0:
        print 'no peaks to subtract from control'
        return

    for observation in plots_in:
        control_ratio = 0.0
        smallest_ratio = 1.0
        for peak_index in plot_to_subtract.peaks:
            height_ratio = observation.y[peak_index] / plot_to_subtract.y[peak_index]
            if height_ratio < smallest_ratio:
                control_ratio = height_ratio
        if control_ratio < 0.0:
            print 'trivial control presence in spectrum, aborting control removal'
            continue
        print 'normalizing control ratio is %s'%(control_ratio)
        new_ys = []
        for i in range(len(observation.y)):
            new_ys.append(observation.y[i] - (control_ratio*plot_to_subtract.y[i]))
        observation.y = new_ys

def standardize_plots(plots_in):
    pass
    """
    # first get any standards information available
    standards = []
    for std_ind in self.project.standards_selected:
        standards.append(string.atof(self.project.standards_weights[std_ind]))
    standards.sort()
    print 'standards', standards
    # so far, coordinates have only been generated, not drawn
    searchrange = .002    # percent weight for tolerance of x- shift
    # right now, this code assumes only one standard is selected
    
    if len(standards) > 1:
        print 'only working with one standard, deleting the others'
        standards = standards[0:0]
    standards_heights = []
    standards_indices = []
    standards_xs      = []
    for mspanel in mspanels:
        standards_heights.append([])       # each panel gets a list of heights for each standard
        standards_indices.append([])       # each panel gets a list of heights for each standard
        standards_xs.append([])
        for stdwt in standards:
            # first locate the x-index of the weight
            midind = -1
            midx   = 0.0
            xs = mspanel.x           # dereference speedup
            for x_ind in range(len(xs)):
                if xs[x_ind] > stdwt:
                    if abs(xs[x_ind]-stdwt) < abs(xs[x_ind-1]-stdwt):
                        midind = x_ind
                        midx   = xs[x_ind]
                    else:
                        midind = x_ind-1
                        midx   = xs[x_ind-1]
                    break
            startind = -1
            startx   = 0.0
            for x_ind in range(len(xs)):
                if midx-xs[x_ind] <= stdwt * searchrange:
                    if abs(xs[x_ind]-stdwt) < abs(xs[x_ind-1]-stdwt):
                        startind = x_ind
                        startx   = xs[x_ind]
                    else:
                        startind = x_ind-1
                        startx   = xs[x_ind-1]
                    break
            endind = -1
            endx   = 0.0
            for x_ind in range(len(xs)):
                if xs[x_ind]-midx > stdwt * searchrange:
                    if abs(xs[x_ind]-stdwt) < abs(xs[x_ind-1]-stdwt):
                        endind = x_ind
                        endx   = xs[x_ind]
                    else:
                        endind = x_ind-1
                        endx   = xs[x_ind-1]
                    break
            maxind = 0
            maxht  = 0.0
            ys = mspanel.y
            for x_ind in range(startind, endind):
                if ys[x_ind] > maxht:
                    maxht = ys[x_ind]
                    maxind = x_ind
            standards_heights[-1].append(string.atof(maxht))
            standards_indices[-1].append(maxind)
            standards_xs[-1].append(mspanel.x[maxind])
    # now find the lowest
    low_height = 100000.0
    for std in standards_heights:
        if std[0] < low_height:
            low_height = std[0]
    # normalize the y-axis
    panel_index = 0
    for mspanel in mspanels:
        ys = mspanel.y
        for y_ind in range(len(ys)):
            ys[y_ind] = ys[y_ind] * (low_height/standards_heights[panel_index][0])
        panel_index += 1
    # and normalize the x-axis
    panel_index = 0
    for mspanel in mspanels:
        dif = mspanel.x[standards_indices[panel_index][0]] - standards[0]            # find the offset of the panel
        xstd = mspanel.x[standards_indices[panel_index][0]]
        for x_ind in range(len(mspanel.x)):
            #mspanel.x[x_ind] = mspanel.x[x_ind] - (mspanel.x[x_ind]/xstd)*dif
            mspanel.x[x_ind] = mspanel.x[x_ind] - (mspanel.x[x_ind]/xstd)*dif
        panel_index += 1
    """


class Plot:
    def __init__(self, source_file=None):
        self.x = []
        self.y = []
        self.source_file = ""
        self.peaks = []          # a list of indices
        self.troughs = []
        self.statistics = {}
        self.label = None
        self.filter_resolution = 0.005
        self.min_peak_height=None
        self.threshold_peak_slope = 0.25
        self.threshold_peak_slope_window_size = 2
        self.mass_range = [0,10000]
        self.peak_style = 'slope'
        if source_file:
            self.load_csv(source_file)
            
    def set_peak_style(self, peak_style):
        self.peak_style = peak_style
        
    def get_peak_style(self):
        return self.peak_style

    def set_mass_range(self, mass_range):
        self.mass_range = mass_range

    def get_mass_range(self):
        return self.mass_range

    def get_filter_resolution(self):
        return self.filter_resolution

    def set_filter_resolution(self, new_resolution):
        self.filter_resolution = new_resolution

    def get_threshold_peak_parms(self):
        """ returns minslope and windowsize for delimiting peaks
        """
        return self.threshold_peak_slope, self.threshold_peak_slope_window_size

    def set_threshold_peak_parms(self, parms):
        """ takes a 2 element parms array with minslope and windowsize
        """
        self.threshold_peak_slope = parms[0]
        self.threshold_peak_slope_window_size = parms[1]
            
    def _recalculate_baseline(self, selection):
        if selection == None or selection == []:
            selection = [[self.x[0], self.x[-1]]]
        sum = 0.0
        count = 0.0
        for cur_range in selection:
            if cur_range[0] > cur_range[1]:
                t_range = cur_range[0]
                cur_range[0] = cur_range[1]
                cur_range[1] = t_range
            for i in range(len(self.x)):
                if self.x[i] > cur_range[0] and self.x[i] < cur_range[1]+1.0:
                    sum += self.y[i]
                    count += 1.0
                    
        if count <= 1:
            print 'no region selected for assessing baseline'
            return
        
        average = sum / count

        accum = 0.0
        for cur_range in selection:
            for i in range(len(self.x)):
                if self.x[i] > cur_range[0] and self.x[i] < cur_range[1]+1.0:
                    accum += (average - self.y[i])**2

        self.peak_SD = math.sqrt(accum/(count-1.0))

        self.peak_LOD =  3.0 * self.peak_SD    # limit of detection
        self.peak_LOQ = 10.0 * self.peak_SD    # limit of quantification
        name = os.path.split(self.source_file)
        if len(name)== 2:
            name = name[1]
        
        print '%s - SD %5.3f, LOD %5.3f, LOQ %5.3f'%(self.source_file, self.peak_SD, self.peak_LOD, self.peak_LOQ)
        self.adjust_min_peak_height(self.peak_LOQ)
        return self.peak_LOQ
    
    def calculate_statistics(self):
        """ calculates .statistics dict, including 'average', 'standard_deviation', 'zscore'.
            z_scores takes up a little memory, might want to delete after use.
        """
        sum = 0
        count = len(self.y)
        for zz in self.y:
            sum = sum + float(zz)
        self.statistics['average'] = sum / count
        avght = self.statistics['average']
        sum = 0
        for zz in self.y:
            sum = sum + pow(float(zz) - avght, 2)
        self.statistics['standard_deviation'] = math.sqrt(sum/(count-1))
        sd = self.statistics['standard_deviation']
        self.statistics['z_scores'] = []
        for zz in self.y:
            self.statistics['z_scores'].append((float(zz)-avght)/sd)
            
    def load_csv_by_GUI(self):
        filename = tkFileDialog.askopenfilename(title = 'Select a project file', defaultextension='.csv', filetypes=[("comma-delimited plot", "*.csv"),("all files", "*")])
        self.load_csv(filename)

    def load_from_array(self, array, array2=None):
        """ load from a nX2 array or 2 nX1 arrays if array2 is present
        """
        tmp_file = open('tmp3946.csv', 'w')
        if array2:
            for i in range(len(array)):
                line = '%s, %s\n'%(array[i], array2[i])
                tmp_file.write(line)
        else:
            for i in range(len(array)):
                line = '%s, %s\n'%(array[i][0], array[i][1])
                tmp_file.write(line)
        tmp_file.close()
        self.load_csv('tmp3946.csv')
        os.remove('tmp3946.csv')
        
    def load_csv(self, csv_file):
        path = os.path.split(csv_file)
        print 'loading plot %s'%(path[-1])
        self.x = []
        self.y = []
        self.source_file = csv_file
        # convert scientific notation to standard
        source = open(csv_file)
        lines = source.readlines()
        for zz in lines:
            if (string.find(zz, "Intensity") > 0):
                continue
            if len(zz) == 0:
                break
            coors = string.split(zz, ',')
            if string.find(coors[0], 'e') > 0:    # check this e stuff for bugs
                self.x.append(float('0.000'))
            else:
                self.x.append(float(coors[0]))
            if string.find(coors[1], 'e') > 0:
                self.y.append(float('0.000'))
            else:
                self.y.append(float(coors[1]))
        source.close()

    def smooth_plot(self, smoothing_gap=None):
        tmpy = 0
        x_temp = self.x
        self.x = []
        y_temp = self.y
        self.y = []
        if not smoothing_gap:
            smoothing_gap = parms.get('smoothing_gap')
        for zz in range(smoothing_gap, len(x_temp)-smoothing_gap-1):
            self.x.append(x_temp[zz])                                 # store the original x-coor
            tmpy = 0.00
            if (smoothing_gap > 0):                                   # skip this loop if no smoothing
                for yy in range(zz-smoothing_gap, zz+smoothing_gap):    # over the range to smooth
                    tmpy = tmpy + y_temp[yy]                # collect the values
                tmpy = tmpy /((smoothing_gap*2)-1)                      # and divide by the width
            else:                                                     # if no smoothing
                tmpy = y_temp[zz]                                       # just copy the y-val directly
            self.y.append(tmpy)
            
    def write(self, filename):
        outfile = open(filename, 'w')
        lines = []
        for i in range(len(self.x)):
            lines.append('%s, %s\n'%(self.x[i], self.y[i]))
        outfile.writelines(lines)
        outfile.close()

    def adjust_min_peak_height(self, min_peak_height):
        self.min_peak_height = min_peak_height
        
    def recognize_peaks(self, min_weight=1000, min_peak_height=None):
        direction = -1          # start out looking for a trough
        last_peak = 0
        last_trough = 0
        self.peaks = []
        self.troughs = []
        self.peak_count = 0
        self.trough_count = 0
        if min_weight == None:
            min_weight = self.get_mass_range()[0]
        height_thresh = parms.get('height_thresh')
        if min_peak_height == None:
            self.min_peak_height = parms.get('min_peak_height')
        else:
            self.min_peak_height = min_peak_height
        # this code performs the recognition of peaks in the sum window
        for i in range(0,len(self.x)):
            if self.x[i] > min_weight:
                if direction == 0:                                                      # if no direction yet
                    if self.y[last_peak] >= self.y[i] + height_thresh:                  #   if last peak >> current
                        direction = -1                                                  #     then decreasing
                    elif self.y[i] >= self.y[last_trough] + height_thresh:              #   else if current >> last trough
                        direction = 1                                                   #     then increasing
                    if self.y[last_peak] < self.y[i]:                                   #   if last peak < current 
                        last_peak = i                                                   #     last peak = current
                    elif self.y[i] < self.y[last_trough]:                               #   else if current < last trough
                        last_trough = i                                                 #     last trough = current
                elif direction == 1:                                                    # else if increasing
                    if self.y[last_peak] < self.y[i]:                                   #   if last peak < current
                        last_peak = i                                                   #     last peak = current
                    elif self.y[last_peak] >= self.y[i] + height_thresh:                #   else if last peak >> current
                        direction = -1                                                  #     direction decreasing
                        last_trough = i                                                 #     last trough = current
                        if self.y[i] > self.min_peak_height:                            #     if current > min peak height
                            self.peaks.append(last_peak)                                #       record this peak
                            self.peak_count = self.peak_count+1                         #
                elif direction == -1:                                                   # else if decreasing
                    if self.y[last_trough] > self.y[i]:                                 #   if last trough > current
                        last_trough = i                                                 #     last trough = current
                    elif self.y[i] >= self.y[last_trough] + height_thresh:              #   else if current >> last trough
                        direction = 1                                                   #     direction increasing
                        last_peak = i                                                   #     last peak = current
                        self.troughs.append(last_trough)                                #     record this trough

        new_troughs = []      # this gets a set of lists, two indices each, representing the bounds of a peak
        for peak in self.peaks:
            # get the next and last troughs
            d = 10000
            last_trough = 0
            for trough in self.troughs:
                dist = abs(self.x[peak] - self.x[trough])
                if dist < d and self.x[trough] < self.x[peak]:
                    d = dist
                    last_trough = trough
            # now search for a close-to-zero point between the peak and the trough
            for j in range(peak, last_trough, -1):
                if self.y[j] < 0.2:
                    last_trough = j
                    break

            d = 10000
            next_trough = 0
            for trough in self.troughs:
                dist = abs(self.x[peak] - self.x[trough])
                if dist < d and self.x[trough] < self.x[peak]:
                    d = dist
                    next_trough = trough
            # now search for a close-to-zero point between the peak and the trough
            for j in range(peak, next_trough):
                if self.y[j] < 0.1:
                    next_trough = j
                    break
            new_troughs.append([last_trough, next_trough])
        self.troughs = new_troughs
        self.trough_count = len(new_troughs)

    def get_peak_count(self):
        return self.peak_count
    
    def optimize_peaks_fmin(self, type):
        self.quantify_peaks(1, 'simplex', type)
        
    def optimize_peaks_powell(self, type):
        self.quantify_peaks(1, 'powell', type)
        
    def approximate_peaks(self, type):
        self.quantify_peaks(0, 0, type)
        
    def fit_peak_approximations(self, fit_type, min_type='simplex'):
        """ assumes recognize peaks was already called """
        self.distribution = Distribution(fit_type)
        peak_args = []
        output = []
        coor_1200 = 0
        for coor_index in range(0,len(self.x)):
            if self.x[coor_index] > self.get_mass_range()[0]:
                coor_1200 = coor_index
                break
        seed = time.clock()
        random.seed(seed)
        #self.distribution.randomize_initial_constants()
        print "optimizing for randomized constants %s"%(self.distribution.initial_constants)
        peak_args = self.distribution.get_parameter_approximations(self)
        print 'initial fitness = %s'%(self.get_distribution_initial_fitness(peak_args, self, self.peaks))
        indices = []
            
        c_lambda = lambda args: self.get_all_panel_fitness(args, coor_1200, len(panel.x))
        if min_type == 'simplex':
            opt_args = scipy.stats.optimize.fmin(c_lambda, self.distribution.initial_constants, args=(), xtol=1e-2, ftol=1e-2, maxiter=1000, full_output=0, disp=1)
        if min_type == 'powell':
            opt_args = scipy.stats.optimize.fmin_powell(c_lambda, self.distribution.initial_constants, args=(), xtol=1e-2, ftol=1e-2, maxiter=100, full_output=0, disp=1)
        if min_type == 'fsolve':
            opt_args = scipy.stats.optimize.fsolve(c_lambda, self.distribution.initial_constants, args=(),fprime=None,full_output=0,col_deriv=0,xtol=1e-2)
        print opt_args, seed

    def get_fitness(self, constants, start, end, verbose=0):
        score = 0.0
        panel_count = 0
        peak_args = self.distribution.get_parameter_approximations(self, constants)
        # first calculate the distributions
        self.distribution.calculate_distributions(peak_args, self, self.peaks, start, end)
        # given this set of parameters, get the sum squared error
        dif = 0.0
        corcnt = 0
        for coor_index in range(start, end):
            if self.x[coor_index] > 3500 and self.x[coor_index] < 4400:
                continue
            else:
                if self.y[coor_index] > 0.5:
                    dif = dif + ((self.y[coor_index]-self.dist_y[coor_index])**2)
                    corcnt += 1
        if verbose:
            print "fitness %11.8f"%(math.sqrt(dif/(0.0+corcnt)))
        score +=  math.sqrt(dif/(0.0+corcnt))
        panel_count = panel_count + 1

        if verbose:
            print "average_fitness = %11.8f"%(score),
            for i in constants:
                print "%5.3f"%(i),
            print
        return score
    
    def quantify_peaks(self, do_optimize, optimization_type=0, fit_type='double gaussian'):
        self.distribution = Distribution(fit_type)
        pre_fitness = []
        post_fitness = []
        fitness_time = []
        # Initialize the scale_y with full length
        peak_args = self.distribution.get_parameter_approximations(self)
                
        # if activated, find a set of parameters for each distribution
        if do_optimize:
            # using these approximations, scan for regions standing out from the baseline
            region_count  = 0
            region_starts = []
            region_ends   = []
            in_region = 0
            index = 0
            region_dividing_thresh = 0.5
            for ht in self.dist_y:
                if not in_region:
                    if ht > region_dividing_thresh:                          # start a region
                        start = index
                        in_region = 1
                else:
                    if ht < region_dividing_thresh:                          # close a region
                        if index-start > 10:                 # if the region is big enough
                            region_starts.append(start)
                            region_ends.append(index)
                            in_region = 0
                            region_count = region_count + 1
                index = index+1
            print "%d regions to optimize"%(region_count)
            portion_sum = 0.0
            for i in range(0,region_count):
                portion = (0.0+region_ends[i] - region_starts[i])/len(self.dist_y)
                print "%d %d - %5.3f pcnt of space"%(region_starts[i], region_ends[i], portion)
                portion_sum = portion_sum + portion
            print "%5.3f pcnt of space covered"%(portion_sum)
    
            for region_index in range(0,region_count):
                print "\nworking on region %d from %6.1f to %6.1f"%(region_index, self.x[region_starts[region_index]], self.x[region_ends[region_index]])
                # collect the peak indices of any peaks in this region
                peak_index = 0
                peak_indices = []
                initial_args = []
                for peak in self.peaks:
                    if peak >= region_starts[region_index] and peak <= region_ends[region_index]:
                        peak_indices.append(peak)
                        print "peak_args (%s)"%(self.peak_args[peak_index])
                        initial_args.append(self.peak_args[peak_index])
                    peak_index = peak_index + 1
                arg_list = []
                for args in initial_args:
                    for arg in args:
                        arg_list.append(float(arg))
                print "contains %s peaks"%(len(peak_indices))
                if len(peak_indices)==0:
                    continue                # to the next region
                # initialize the distributions
                self.get_distribution_initial_fitness(initial_args, self, peak_indices)
                # create a fitness function for this region
                fit_func = lambda args: self.get_extended_distribution_fitness(args, self, peak_indices, region_starts[region_index], region_ends[region_index])
                # perform the optimization
                fmin_max_iter = 500
                powell_max_iter = 10
                anneal_max_iter = 100
                if optimization_type == 'simplex':
                    opt_args = scipy.stats.optimize.fmin(fit_func, arg_list, args=(), xtol=1e-2, ftol=1e-2, maxiter=fmin_max_iter, full_output=0, disp=1)
                elif optimization_type == 'powell':
                    opt_args = scipy.stats.optimize.fmin_powell(fit_func, arg_list, args=(), xtol=1e-2, ftol=1e-2, maxiter=powell_max_iter, full_output=0, disp=1)
                #elif optimization_type == 2:
                #    opt_args = scipy.stats.optimize.anneal(panel_fit_func, initial_args, upper=1000.0, lower=0.0, feps=1e-4, maxiter=anneal_max_iter, schedule='fast')
                #post_fitness.append(self.get_extended_distribution_fitness(opt_args, panel, peak_indices, region_starts[region_index], region_ends[region_index]))
                print "args %s"%(opt_args)
                for cnt in range(0,len(peak_indices)):                              # go through each of the peaks included in this region
                    for peak_index in range(0,len(self.peaks)):                    # search through the peaks
                        if self.peaks[peak_index] == peak_indices[cnt]:            # if this is the one
                            for i in range(self.distribution.arg_cnt):
                                self.peak_args[peak_index][i] = opt_args[(cnt*self.distribution.arg_cnt)+i]         # set the SD for this peak to the SD from the optimization output
        self.distribution.plot_distributions(self)

    def get_approximations_fitness(self, input_args, start, end):
        args = self.distribution.get_parameter_approximations(self, input_args,fit_type)
        soln = self.get_distribution_fitness(args, self, self.peaks, start, end)
        return soln
    
    def get_extended_distribution_fitness(self, peak_args, peak_indices, start, end):
        # turns a flat list of args into a list of lists
        long_list = []
        cnt = 0
        for arg in range(0,len(peak_indices)):
            short_list = []
            for i in range(0,self.distribution.arg_cnt):
                short_list.append(peak_args[cnt])
                cnt = cnt + 1
            long_list.append(short_list)
        return self.get_distribution_fitness(long_list, self, peak_indices, start, end)

    def get_distribution_initial_fitness(self, peak_args, panel, peak_indices):
        self.distribution.calculate_initial_distributions(peak_args, self, peak_indices)
        start = 0
        end = len(self.x)
        dif = 0.0
        for coor_index in range(start, end):
            dif = dif + ((self.y[coor_index]-self.dist_y[coor_index])**2)
        return math.sqrt(dif/(end-start-1.0))
        
    def get_distribution_fitness(self, peak_args, panel, peak_indices, start, end):
        # first calculate the distributions
        self.distribution.calculate_distributions(peak_args, self, peak_indices, start, end)
        # given this set of parameters, get the sum squared error
        dif = 0.0
        for coor_index in range(start, end):
            dif = dif + ((self.y[coor_index]-self.dist_y[coor_index])**2)
        return math.sqrt(dif/(end-start-1.0))
    def set_label(self, label):
        # this should be a float
        self.label = label
    def get_label(self):
        return self.label
    def hard_rescale_y(self, ratio):
        new_ys = []
        ratio = ratio/100.0
        for yval in self.y:
            new_ys.append(ratio * yval)
        self.y = new_ys

    def get_area_for_peak_index(self, peak, style='slope'):
        # actually just calculates the sum of peak heights, but this should
        # be proportional to the area, since the units don't really change
        # in local regions
        # this can be modified by simply dividing by the x-range
        a,b,c = self.get_trough_indices_for_peak_index(peak, style)
        sum = 0.0
        for i in range(a,c+1):
            sum += self.y[i]
        return sum
            
    def get_trough_indices_for_peak_index(self, peak, style='slope'):
        """
        # search for the highest point within the resolution of the machine        
        max_height = 0.0
        tallest = None
        for i in range(len(self.y)):
            d = abs(self.x[i] - self.x[peak])
            if d < self.filter_resolution * self.x[peak]:
                if self.y[i] > max_height:
                    max_height = self.y[i]
                    tallest = i
        if tallest == None:
            closest = peak
        else:
            closest = tallest
        """
        # search for the closest peak within the resolution of the machine
        closest_peak_distance = 100000
        closest_peak = None
        for peak_index in self.peaks:
            d = abs(self.x[peak_index]-self.x[peak])
            if d < closest_peak_distance:
                if d < self.filter_resolution * self.x[peak]:
                    closest_peak_distance = d
                    closest_peak = peak_index
        if closest_peak == None:
            closest = peak
        else:
            closest = closest_peak

        if self.get_peak_style() == 'minrange':
            # get the last trough from those originally defined
            idrange = 25.0   # a maximum of idrange from the peak
            d = 10000
            last_trough = 0
            size_thresh = 0.2
            for trough in self.troughs:  # search for the closest
                dist = abs(self.x[closest] - self.x[trough])
                if dist < idrange and self.x[trough] < self.x[closest]:
                    d = dist
                    last_trough = trough
            if last_trough == 0:
                # none was assigned. Search back for the closest index to 5.0
                for i in range(1000):
                    if self.x[closest-i] > self.x[closest] - idrange:
                        last_trough = closest-i
                    else:
                        break
    
            # now search for a close-to-zero point between the peak and the trough
            for j in range(closest, last_trough, -1):
                if self.y[j] < size_thresh:
                    last_trough = j
                    break

            d = 10000
            next_trough = 0
            for trough in self.troughs:
                dist = abs(self.x[closest] - self.x[trough])
                if dist < idrange and self.x[trough] > self.x[closest]:
                    d = dist
                    next_trough = trough
    
            if next_trough == 0:
                # none was assigned. Search forward for the closest index to 2*peakheight
                for i in range(1000):
                    if self.x[closest+i] < self.x[closest] + idrange:
                        next_trough = closest+i
                    else:
                        break
    
            # now search for a close-to-zero point between the peak and the trough
            for j in range(closest, next_trough):
                if self.y[j] < size_thresh:
                    next_trough = j
                    break
                
        elif self.get_peak_style() == 'slope':
            minslope    = self.threshold_peak_slope      # must be below a minimum slope
            windowsize  = self.threshold_peak_slope_window_size
            last_trough = closest
            next_trough = closest
            min_height  = 0.75     # must be below a certain fraction of the peak height
            searching_for_downslope = 1
            for i in range(0, 100):
                if searching_for_downslope and self.y[closest-i] > self.y[closest-i+1]:
                    continue
                searching_for_downslope = 0
                if self.y[closest-i] < min_height*self.y[closest]:
                    # find the average slope over the next windowsize x-values
                    slope = 0.0
                    for j in range(0,windowsize):
                        slope += (self.y[closest-i-j]-self.y[closest-i-j-1]) / (self.x[closest-i-j]-self.x[closest-i-j-1])
                    slope /= float(windowsize)
                    if slope > minslope:
                        last_trough = closest-i
                    else:
                        break
                else:
                    firstslope = (self.y[closest-i]-self.y[closest-(i+1)]) / (self.x[closest-i]-self.x[closest-i-1])
                    if firstslope < 0.0:
                        last_trough = closest-i
                        break

            searching_for_downslope = 1
            for i in range(100):
                if searching_for_downslope and self.y[closest+i] > closest+i+1:
                    continue
                searching_for_downslope = 0
                if self.y[closest+i] < min_height*self.y[closest]:
                    # find the average slope over the next windowsize x-values
                    slope = 0.0
                    for j in range(0,windowsize):
                        slope += (self.y[closest+i+j]-self.y[closest+i+j+1]) / (self.x[closest+i+j+1]-self.x[closest+i+j])
                    slope /= float(windowsize)
                    if slope > minslope:
                        next_trough = closest+i
                    else:
                        break
                else:
                    firstslope = (self.y[closest+i]-self.y[closest+i+1]) / (self.x[closest+i+1]-self.x[closest+i])
                    if firstslope < 0.0:
                        next_trough = closest+i
                        break

        #if last_trough == closest or next_trough == closest:
        #    print '### unshifted bounds b%5.3f m%5.3f e%5.3f'%(self.x[last_trough], self.x[closest], self.x[next_trough])
        if last_trough == closest:
            last_trough = closest-1
        if next_trough == closest:
            next_trough = closest+1
        return last_trough, closest, next_trough
    
