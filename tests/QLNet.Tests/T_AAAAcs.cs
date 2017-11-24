/*
 Copyright (C) 2008 Siarhei Novik (snovik@gmail.com)
 Copyright (C) 2008 Andrea Maggiulli

 This file is part of QLNet Project https://github.com/amaggiulli/qlnet

 QLNet is free software: you can redistribute it and/or modify it
 under the terms of the QLNet license.  You should have received a
 copy of the license along with this program; if not, license is
 available online at <https://github.com/amaggiulli/qlnetLicense.html>.

 QLNet is a based on QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/
 The QuantLib license is available online at http://quantlib.org/license.shtml.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

using System;
using System.Collections.Generic;
using System.Linq;
#if NET40 || NET45
using Microsoft.VisualStudio.TestTools.UnitTesting;
#else
   using Xunit;
#endif
using QLNet;

namespace TestSuite
{
#if NET40 || NET45
   [TestClass()]
#endif
   class T_AAAAcs
    {
#if NET40 || NET45
      [TestMethod()]
#else
       [Fact]
#endif
      public void TestMethod1()
      {
         double a = 3;
         Assert.IsTrue(a == 3);
      }

      [TestMethod]
      public void TestMethod2()
      {
         Assert.IsTrue(false, "hello");
      }
   }
}
