import os
from time import sleep

while(os.system("matlab -nodisplay")):
   print "waiting for license..."
   sleep(1)
