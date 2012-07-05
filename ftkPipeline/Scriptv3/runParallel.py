#!/usr/bin/env python
# -*- coding: utf-8 -*-

import threading
import shutil
import fnmatch
import os
import subprocess
import os.path
import time
import glob
import run

from multiprocessing import Process


DATA_FOLDER_ALL = ['/0113_NRRD_CROPPED','/0117_NRRD_CROPPED','/0120_NRRD_CROPPED','/0123_NRRD_CROPPED','/0128_NRRD_CROPPED','/0131_NRRD_CROPPED','/0323_NRRD_CROPPED','/0405_NRRD_CROPPED','/0409_NRRD_CROPPED','/0410_NRRD_CROPPED','/0412_NRRD_CROPPED','/1206_NRRD_CROPPED']

#DATA_FOLDER_ALL = ['/0131_test','/0131_test2']


#p1 = Process(target=run.main, args=('/0113_NRRD_CROPPED',))
#p2 = Process(target=run.main, args=('/0113_NRRD_CROPPED',))
#p1 = Process(target=run.main, args=('/0131_test',))
#p2 = Process(target=run.main, args=('/0131_test2',))
#p1.start()
#p2.start()
#p.join()

for DATA_FOLDER in DATA_FOLDER_ALL:
	p1 = Process(target=run.main, args=(DATA_FOLDER,))
	p1.start()

#for DATA_FOLDER in DATA_FOLDER_ALL:
	#threading.Thread(target=run.main(DATA_FOLDER)).start()
