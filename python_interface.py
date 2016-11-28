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
#Modify the below line for your own machine
pathToPrevariety = "/home/jeff/Desktop/HybridPrevariety/"

def TropicalPrevariety(support):
	support = str(support)
	support = support.replace("], ", "]")
	support = support.replace(" ","")
	prevariety = Popen(
		pathToPrevariety + "prevariety.out", stdin=PIPE, stdout=PIPE, stderr=PIPE, shell=true)
	ans, err = prevariety.communicate(input = support)
	if len(ans) > 0 and ans[0] != '[':
		raise Exception("Internal error in tropical_prevariety")
	return eval(ans)
