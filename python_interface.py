"""
	#Cyclic 4
	R.<x1,x2,x3,x4> = QQ[]
	polys = [x1+x2+x3+x4,x1*x2+x2*x3+x3*x4+x4*x1,x1*x2*x3+x2*x3*x4+x3*x4*x1+x4*x1*x2]
	support = [[[Integer(j) for j in i] for i in poly.exponents()] for poly in polys]
	TropicalPrevariety(support)
	#Should be equivalent (up to homogenization) to:
	R.ideal(polys).groebner_fan().tropical_intersection().rays()

	#Reduced cyclic 8
	R.<y_1,y_2,y_3,y_4,y_5,y_6,y_7> = QQ[]
	polys = [1 + y_1 + y_2 + y_3 + y_4 + y_5 + y_6 + y_7,y_1 + y_1*y_2 + y_2*y_3
	+ y_3*y_4 + y_4*y_5 + y_5*y_6 + y_6*y_7 + y_7,y_1*y_2 + y_1*y_2*y_3
	+ y_2*y_3*y_4 + y_3*y_4*y_5 + y_4*y_5*y_6 + y_5*y_6*y_7
	+ y_6*y_7 + y_7*y_1,y_1*y_2*y_3 + y_1*y_2*y_3*y_4 + y_2*y_3*y_4*y_5 
	+ y_3*y_4*y_5*y_6 + y_4*y_5*y_6*y_7 + y_5*y_6*y_7 + y_6*y_7*y_1
	+ y_7*y_1*y_2,y_1*y_2*y_3*y_4 + y_1*y_2*y_3*y_4*y_5 + y_2*y_3*y_4*y_5*y_6
	+ y_3*y_4*y_5*y_6*y_7 + y_4*y_5*y_6*y_7 + y_5*y_6*y_7*y_1 + y_6*y_7*y_1*y_2
	+ y_7*y_1*y_2*y_3,y_1*y_2*y_3*y_4*y_5 + y_1*y_2*y_3*y_4*y_5*y_6
	+ y_2*y_3*y_4*y_5*y_6*y_7 + y_3*y_4*y_5*y_6*y_7 + y_4*y_5*y_6*y_7*y_1
	+ y_5*y_6*y_7*y_1*y_2 + y_6*y_7*y_1*y_2*y_3
	+ y_7*y_1*y_2*y_3*y_4,y_1*y_2*y_3*y_4*y_5*y_6 + y_1*y_2*y_3*y_4*y_5*y_6*y_7
	+ y_2*y_3*y_4*y_5*y_6*y_7+ y_3*y_4*y_5*y_6*y_7*y_1 + y_4*y_5*y_6*y_7*y_1*y_2
	+ y_5*y_6*y_7*y_1*y_2*y_3+ y_6*y_7*y_1*y_2*y_3*y_4 + y_7*y_1*y_2*y_3*y_4*y_5]
	support = [[[Integer(j) for j in i] for i in poly.exponents()] for poly in polys]
	TropicalPrevariety(support)
"""

from subprocess import Popen, PIPE
#The below should work for generic machines
import os, inspect
pathToPrevariety = os.path.dirname(inspect.stack()[0][1]) + '/'

def TropicalPrevariety(support, ProcessCount = 1):
	support = str(support)
	support = support.replace("], ", "]")
	support = support.replace(" ","")
	if ProcessCount < 10:
	   support = '0' + str(ProcessCount) + support
	else: 
	   support = str(ProcessCount) + support
	
	prevariety = Popen(
		pathToPrevariety + "prevariety", stdin=PIPE, stdout=PIPE, stderr=PIPE, shell=true)
	ans, err = prevariety.communicate(input = support)
	print(ans)
	#if len(ans) > 0 and ans[0] != '[':
	#	raise Exception("Internal error in tropical_prevariety")
	return
