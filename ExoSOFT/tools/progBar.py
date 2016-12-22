from __future__ import absolute_import
import sys

class ProgBar(object):
    """
    Call in a loop to create terminal progress bar
    
    Heavily modified, but code originally copied from:
    http://stackoverflow.com/questions/3173320/text-progress-bar-in-the-console
    
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        barLength   - Optional  : character length of bar (Int)
    """
    def __init__(self,total=100, decimals = 0, barLength = 100):
        self.formatStr = "{0:." + str(decimals) + "f}"
        self.total = total
        self.decimals = decimals
        self.barLength = barLength
        
    def render(self,iteration,prefix = '', suffix = ''):
        percents        = self.formatStr.format(100 * (iteration / float(self.total)))
        filledLength    = int(round(self.barLength * iteration / float(self.total)))
        bar=''
        if self.barLength>0:
            if filledLength<2:
                filledLength=1
            bar = '=' * (filledLength-1)+'>' + '-' * (self.barLength - filledLength)
        sys.stdout.write('\r%s %s %s%s %s' % (prefix, bar, percents, '%', suffix))
        sys.stdout.flush()
        if iteration == self.total:
            #sys.stdout.write('\n')  #$$  Might want to un-comment this in the future
            sys.stdout.flush()