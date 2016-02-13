#!/usr/bin/env python
import unittest
import math
import orbit
import timeit
import numpy as np

# extra function for speed testing
def speedTestPy(ntimes):
    for i in range(0,ntimes):	
        a = math.sin(1.0)

class TestSwig(unittest.TestCase):
    """
    Test linking with swig. Use simple cases to mimic the functions 
    we will encounter in this package. If this test fail, then the 
    package is not portable on the local machine. 
    """
#     #tests for hello.cc
#     def testInt(self):
#         #pass a int from python as an argument of a c++ function, return an int
#         self.assertEqual(hello.fact(4),24)
#     def testFloat(self):
#         #pass a float from python as an argument of a c++ function, return an
#         #float.
#         self.assertAlmostEqual(hello.csin(0.5*math.pi),math.sin(0.5*math.pi))
#     def testStr(self):
#         #pass a string from python as an argument of a c++ function.
#         hello.speak("helloworld")
#         self.assertAlmostEqual(1,1)
#     def testbool(self):
#         #pass a bool from python as an argument of a c++ function.
#         a = hello.selection(True)
#         self.assertEqual(a,1)
#         b = hello.selection(False)
#         self.assertEqual(b,0)
#     def testNparray1d_scalar(self):
#         #pass an 1-d numpy array from python as an argument of a c++ function.
#         #the value of the array is not changed in the c++ function.
#         a = np.arange(50)
#         self.assertAlmostEqual(hello.sum_scalar(a,1),np.sum(a)+50)
#     def testNparray1d_1d(self):
#         #pass two 1-d numpy array from python as an argument of a c++ function.
#         #the value of the array is not changed in the c++ function.
#         a = np.arange(50)
#         self.assertAlmostEqual(hello.sum_1d(a,a),np.sum(a)*2)
#     def testNparray2d(self):
#         #pass a 2-d numpy array from python as an argument of a c++ function.
#         #the value of the array is not changed in the c++ function.
#         a = np.array([[1,2],[3,4]])
#         self.assertAlmostEqual(hello.sum2d(a),10)
#     # tests for initial version of orbCalc.cpp
#     def testNparray2d(self):
#         #pass a 2-d numpy array from python as an argument of a c++ function.
#         #the value of the array is not changed in the c++ function.
#         print '\n'
#         a = np.array([[1.0,2.0],[3.0,4.0]])
#         #print "\na.shape = "+repr(a.shape)
#         self.assertAlmostEqual(orbCalc.sum2d_double(a),10.0)
#     def testModDouble(self):
#         a = 10.0
#         self.assertAlmostEqual(orbCalc.modAndReturnDouble(a),20.0)
#     def testSpeedTest(self):
#         tic=timeit.default_timer()
#         speedTestPy(1000000)
#         toc=timeit.default_timer()
#         print "\nfor python it took: "+str(toc-tic)
#         tic=timeit.default_timer()
#         orbCalc.speedTestSWIG(1000000)
#         toc=timeit.default_timer()
#         print "for SWIG it took: "+str(toc-tic)
#         self.assertAlmostEqual(10,10)
#     def testNparray2d_modInPlace(self):
#         #pass a 2-d numpy array from python as an argument of a c++ function.
#         #the value of the array is not changed in the c++ function.
#         a = np.zeros((2,2))
#         #print "\na.shape = "+repr(a.shape)
#         #print "a = "+repr(a)
#         orbCalc.modAndReturn2d_double(a)
#         #print "\na.shape = "+repr(a.shape)
#         #print "a = "+repr(a)
#         #both of these tests syntaxs work!
#         np.testing.assert_array_almost_equal(a,np.array([[1.0,2.0],[2.0,3.0]]))
#         #self.assertTrue(np.allclose(a,np.array([[1.0,2.0],[2.0,3.0]])))
#         #print "\na.shape = "+repr(a.shape)
#         #print "a = "+repr(a)
    # tests for initial version of orbit.cc
    def testCPPobj(self):
        print '\nstarting c++ obj test'
        obj = orbit.Orbit()
        obj.testDouble = 20.0
        a = np.array([[1.0,2.0],[3.0,4.0]])
        obj.loadRealData(a)
        obj.loadConstants(10, 20, 30, 40,50,60)
        b = np.zeros((2,2))
        params = np.array([1.0,2.0,3.0,4.0,5.0,6.0])        
        print 'before b = '+repr(b)
        obj.calculate(b,params)
        print 'b = '+repr(b)
        print 'params = '+repr(params)
        np.testing.assert_array_almost_equal(b,np.array([[1.0,2.0],[2.0,3.0]]))

if __name__ == "__main__":
    unittest.main()

