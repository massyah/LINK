import sys
import threading 

class AWorker(threading.Thread):
	def __init__(self,uid):
		print "Will work on file",uid
		self.uid=uid
		threading.Thread.__init__(self)

	def file(self):
		return self.uid

	def run(self):
		my_list=[]
		# Look for prime number from 1 to 1e6
		print "Thread",self.uid,"started working"
		f=open("temp_%d.dat"%(self.uid),"w")
		for n in xrange(1,1000000,1):
			is_prime=True
			for adiv in xrange(2,n,1):
				if n % adiv == 0:
					is_prime=False
					break
			if is_prime:
				my_list.append(n)
				f.write("%d\n"%(n))
		print "Thread",self.uid,"finished working"
		f.close()
		sys.stdout.flush()


threadLock = threading.Lock()
threads = []

# Create new threads
thread1 = AWorker(1)
thread2 = AWorker(2)

# Start them 
thread1.start()
thread2.start()

# Add threads to thread list
threads.append(thread1)
threads.append(thread2)

# Wait for all threads to complete
for t in threads:
    t.join()
print "Exiting Main Thread"

