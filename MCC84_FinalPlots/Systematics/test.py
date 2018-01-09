f1 = open('event_info_sel2_cv.txt','r')
overlap = 0

for l1 in f1.readlines():

  subrun1 = l1.split(' ')[0]
  event1  = l1.split(' ')[1]

  print "Checking event overlap " , overlap

  f2 = open('event_info_sel2_enhancedexttpc.txt','r')
  for l2 in f2.readlines():

    subrun2 = l2.split(' ')[0]
    event2  = l2.split(' ')[1]

    if int(subrun2) == int(subrun1) and int(event2) == int(event1) :
      overlap += 1
      #print "Found one!"
      break

print "Final sample overlap: ", overlap

  
