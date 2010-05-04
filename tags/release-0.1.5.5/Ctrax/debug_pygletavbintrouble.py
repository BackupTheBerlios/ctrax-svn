import ctypes
avbin = ctypes.cdll.avbin
import pyglet.media as media
filename = 'movie20071009_163231_frames0001to0100_xvid.avi'
source = media.load(filename)
ts0 = source.get_next_video_timestamp()
im = source.get_next_video_frame()
seektime = 1./1000000.
source._seek(seektime)
ts0_1 = source.get_next_video_timestamp()
