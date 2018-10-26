#!/usr/bin/env python                                                                                                              
import numpy as np
import conv

run = 'storm/test7'
              
conv.WPE_conv(run)
conv.WPE_conv_direct(run)
conv.WKE_conv(run)
conv.WKE_conv_direct(run)

