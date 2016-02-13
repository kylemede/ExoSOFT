# -*- coding: utf-8 -*-
# Copyright: 2009 Nadia Alramli
# License: BSD
"""Draws an animated terminal progress bar
Usage:
    p = ProgressBar("blue")
    p.render(percentage, message)
"""

import terminal
import sys

class ProgressBar(object):
    """Terminal progress bar class"""
    TEMPLATE = (
     '%(percent)-3s%% %(color)s%(progress)s%(lastblock)s%(normal)s%(empty)s %(message)s\n'
    )
    PADDING = 7

    def __init__(self, color=None, width=None, block='█',
                 lastblock='█', empty=' '):
        """
        color -- color name (BLUE GREEN CYAN RED MAGENTA YELLOW WHITE BLACK)
        width -- bar width (optional)
        block -- progress display character (default '█')
        empty -- bar display character (default ' ')
        """
        if color:
            self.color = getattr(terminal, color.upper())
        else:
            self.color = ''
        if width and width < terminal.COLUMNS - self.PADDING:
            self.width = width
        else:
            # Adjust to the width of the terminal
            self.width = terminal.COLUMNS - self.PADDING
        self.block = block
        self.lastblock = lastblock
        self.empty = empty
        self.progress = None
        self.lines = 0

    def render(self, percent, message=''):
        """Print the progress bar
        percent -- the progress percentage %
        message -- message string (optional)
        """

        # Number of progress bars:
        try:
            percent = list(percent)
            message = list(message)
            n = len(percent)
            if len(message) < n:
                # Pad with blanks
                message += '' * (n - len(message))
        except:
            try:
                percent = [int(percent)]
                message = [message]
                n = 1
            except:
                print "Cannot make a list of inputs to MultiProgressBar.render"
                sys.exit()

        for i in range(n):
            msg = message[i]
            pct = int(percent[i])
            inline_msg_len = 0
            if msg:
                # The length of the first line in the message
                inline_msg_len = len(msg.splitlines()[0])
            if inline_msg_len + self.width + self.PADDING > terminal.COLUMNS:
                # The message is too long to fit in one line.
                # Adjust the bar width to fit.
                bar_width = terminal.COLUMNS - inline_msg_len -self.PADDING
            else:
                if pct > 99:
                    bar_width = self.width
                else:
                    bar_width = self.width + 1

            # Check if render is called for the first time
            if i == 0 and self.progress is not None:
                self.clear()
                self.lines = 0
            self.progress = (bar_width * pct) / 100
            showlast = 0
            if self.progress > 0:
                showlast = 1
                
            data = self.TEMPLATE % {
                'percent': pct,
                'color': self.color,
                'progress': self.block * (self.progress),
                'lastblock': self.lastblock * (showlast),
                'normal': terminal.NORMAL,
                'empty': self.empty * (bar_width - self.progress - showlast),
                'message': msg
                }
            sys.stdout.write(data)
            sys.stdout.flush()
            # The number of lines printed
            self.lines += len(data.splitlines())

    def clear(self):
        """Clear all printed lines"""
        sys.stdout.write(
            self.lines * (terminal.UP + terminal.BOL + terminal.CLEAR_EOL)
        )

