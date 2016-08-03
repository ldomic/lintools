#!/usr/bin/env python
# -*- coding: utf-8 -*- 

from termcolor import colored
import sys

PROGRESS_BLOCKS = u" ▏▎▍▌▋▊█"
PROGRESS_SPINNER = [u"🌘",u"🌗",u"🌖",u"🌕",u"🌔",u"🌓",u"🌒",u"🌑"]

class ProgressBar(object):
    def __init__(self, title, width=30):
        self.title = title
        self.spinner_state = 0
        self.width = width
    def update(self, progress=None, colour="green"):
        """"""
        percentage = progress/1 * 100
        self.spinner_state += 1
        self.spinner_state %= len(PROGRESS_SPINNER)
        spinner_char = PROGRESS_SPINNER[self.spinner_state]
        progress_bar = self.width * progress
        full_blocks = int(progress_bar)
        partial_block = int((progress_bar - full_blocks) * 8)
        blank_blocks = self.width - (full_blocks + 1)
        block_chars = [0 for x in range(self.width)]
        for x in range(full_blocks):
            if x < len(block_chars):
                block_chars[x] = -1
        if full_blocks < self.width:
            block_chars[full_blocks] = partial_block
        progress_bar_text = u"┣" +"".join([PROGRESS_BLOCKS[x] for x in block_chars])+ u"┫"
        percent_text = str(int(percentage))+"%"
        sys.stdout.write(colored(u"   "+spinner_char+"  " + self.title.ljust(30) + progress_bar_text +"   "+percent_text+ "\r" , colour, attrs=['bold']))
        sys.stdout.flush()
    def finish(self, colour="green"):
        self.update(1, colour)
        sys.stdout.write(colored("   ✔  " + self.title.ljust(30) +  "\n", "green", attrs=['bold']))