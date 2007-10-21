import os
from time import sleep

while(os.system("matlab")):
   print "waiting for license..."
   sleep(1)
