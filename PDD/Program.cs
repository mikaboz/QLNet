using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using QLNet;
namespace PDD
{
 
   class Program
   {

      static void MCBarrierngineTest()
      {
         #region General

         Calendar calendar = new TARGET();
         DayCounter dayCounter = new Actual365Fixed();
         Date today = DateTime.Now;
         today = calendar.adjust(today);
         #endregion

         #region BlackScholesProcess

         double spot = 100;
         double flatRate = 0.015;
         double flatDividend = 0;
         double flatVol = 0.45;

         IDiscretization1D discretization = new EulerDiscretization();
         Handle<BlackVolTermStructure> volatilityTS = new Handle<BlackVolTermStructure>(new BlackConstantVol(today, calendar, flatVol, dayCounter));
         Handle<YieldTermStructure> riskFreeTS = new Handle<YieldTermStructure>(new FlatForward(today, flatRate, dayCounter));
         Handle<YieldTermStructure> dividendTS = new Handle<YieldTermStructure>(new FlatForward(today, flatDividend, dayCounter));
         Handle<Quote> hspot = new Handle<Quote>(new SimpleQuote(spot));
         GeneralizedBlackScholesProcess process = new BlackScholesMertonProcess(hspot, dividendTS, riskFreeTS, volatilityTS, discretization);
         #endregion

         #region PDD

         double barrier = 150;
         
         Date maturityDate = calendar.advance(today, 2, TimeUnit.Months);

         BarrierOption barrierOption = new BarrierOption(Barrier.Type.UpOut, barrier, 0, new PlainVanillaPayoff(Option.Type.Call, 110),new EuropeanExercise(maturityDate));
         #endregion

         #region MonteCarlo
         int requiredSamples = 1002;
         ulong seed = 32;
         MCBarriereEngine<PseudoRandom, Statistics> mcBarrierEngine = new MCBarriereEngine<PseudoRandom, Statistics>(process, 100, null, false, true, false, requiredSamples, null, null, seed);
         barrierOption.setPricingEngine(mcBarrierEngine);
         double value = barrierOption.NPV();
         double error = barrierOption.errorEstimate();
         Console.WriteLine(mcBarrierEngine.sampleAccumulator().samples());
         Console.WriteLine(value);
         Console.WriteLine(error);
         Console.ReadLine();
         #endregion
      }

      static void PDDMCEngineTest()
      {
         #region General

         Calendar calendar = new TARGET();
         DayCounter dayCounter = new Actual365Fixed();
         Date today = DateTime.Now;
         today = calendar.adjust(today);
         #endregion

         #region BlackScholesProcess

         double spot = 100;
         double flatRate = 0.015;
         double flatDividend = 0;
         double flatVol = 0.45;

         IDiscretization1D discretization = new EulerDiscretization();
         Handle<BlackVolTermStructure> volatilityTS = new Handle<BlackVolTermStructure>(new BlackConstantVol(today, calendar, flatVol, dayCounter));
         Handle<YieldTermStructure> riskFreeTS = new Handle<YieldTermStructure>(new FlatForward(today, flatRate, dayCounter));
         Handle<YieldTermStructure> dividendTS = new Handle<YieldTermStructure>(new FlatForward(today, flatDividend, dayCounter));
         Handle<Quote> hspot = new Handle<Quote>(new SimpleQuote(spot));
         GeneralizedBlackScholesProcess process = new BlackScholesMertonProcess(hspot, dividendTS, riskFreeTS, volatilityTS, discretization);
         #endregion

         #region PDD

         double acquisition = 120;
         double percent = 0.8;

         Period provFrequency = new Period(Frequency.Annual);
         Date maturityDate = calendar.advance(today, 2,TimeUnit.Months);//provFrequency);
         Period conditionalPeriod = new Period(3, TimeUnit.Months);

         PDD pdd = new PDD(acquisition, percent, maturityDate, conditionalPeriod);
         #endregion

         #region MonteCarlo
         int timeStepsPerMaximumConditionPeriod = 3;
         int requiredSamples = 1002;
         ulong seed = 32;
         MCPDDEngine<PseudoRandom, Statistics> mCPDDEngine =
            new MCPDDEngine<PseudoRandom, Statistics>(process, timeStepsPerMaximumConditionPeriod, false, false, false,requiredSamples,null,null, seed);
         pdd.setPricingEngine(mCPDDEngine);
         double value = pdd.NPV();
         double error = pdd.errorEstimate();
         Console.WriteLine(mCPDDEngine.sampleAccumulator().samples());
         Console.WriteLine(value);
         Console.WriteLine(error);
         Console.ReadLine();
         #endregion
      }

      static void Main(string[] args)
      {
         MCBarrierngineTest();
      }
   }
}
