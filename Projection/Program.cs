using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Projection
{
   using QLNet;
   class Program
   {
      static void Main(string[] args)
      {
         ProjectionBlackScholes();
      }

      static void ProjectionBlackScholes()
      {
         DateTime now = DateTime.Now;

         #region TimeGrid : horizon 10 ans, fréquence journalière
         Calendar calendar = new TARGET();
         Date today = DateTime.Today;
         Period horizon = new Period(10, TimeUnit.Years);
         Date finalDate = today + horizon;
         today = calendar.adjust(today);
         finalDate = calendar.adjust(finalDate,BusinessDayConvention.Preceding);
         DayCounter dayCounter = new Actual365Fixed();
         double timeHorizon = dayCounter.yearFraction(today, finalDate);
         int steps = 100;
         TimeGrid timeGrid = new TimeGrid(timeHorizon, steps);
         Console.WriteLine("timeGrid:" +
            "\nHorizon: " + timeHorizon +
            "\nSteps: " + steps +
            "\nStart:" + timeGrid.First()+
            "\nSize: " + timeGrid.size()
            );
         Console.ReadLine();
         #endregion

         #region BlackScholesProcess
         IDiscretization1D discretization = new EulerDiscretization();
         IParametricVolatility parametricVolatility = new ExponentialVolatilityParameter(0.1, 1, 0.2, 1, horizon);
         Handle<BlackVolTermStructure> volatilityTS = new Handle<BlackVolTermStructure>(new ParametricVolatilityTermStructure(parametricVolatility, today, dayCounter, calendar));
         Handle<YieldTermStructure> riskFreeTS = new Handle<YieldTermStructure>(new FlatForward(today, 0.01, dayCounter));
         Handle<YieldTermStructure> dividendTS = new Handle<YieldTermStructure>(new FlatForward(today, 0.0, dayCounter));
         Handle<Quote> spot = new Handle<Quote>(new SimpleQuote(100));
         StochasticProcess process = new GeneralizedBlackScholesProcess(spot, dividendTS, riskFreeTS, volatilityTS, null,discretization);
         #endregion

         #region RandomSequenceGenerator: MersennTwister + InverseCumulativeNormal
         RandomSequenceGenerator<MersenneTwisterUniformRng>
            URNG = new RandomSequenceGenerator<MersenneTwisterUniformRng>(steps,0);
         InverseCumulativeRsg<RandomSequenceGenerator<MersenneTwisterUniformRng>, InverseCumulativeNormal>
            RSG = new InverseCumulativeRsg<RandomSequenceGenerator<MersenneTwisterUniformRng>, InverseCumulativeNormal>(
               URNG,
               new InverseCumulativeNormal());
         #endregion

         #region Générateur de trajectoires
         PathGenerator<InverseCumulativeRsg<RandomSequenceGenerator<MersenneTwisterUniformRng>, InverseCumulativeNormal>>
            pathGenerator = new PathGenerator<InverseCumulativeRsg<RandomSequenceGenerator<MersenneTwisterUniformRng>, InverseCumulativeNormal>>(
               process,
               timeGrid,
               RSG,
               false);
         #endregion

         int nSimulations = 500;
         SequenceStatistics stats = new SequenceStatistics(steps+1);
         Statistics optionStat = new Statistics();
         for (int i = 0; i < nSimulations; i++)
         {
            Path path = (Path)pathGenerator.next().value;
            List<double> values = new List<double>();
            for (int j = 0; j < path.length(); j++)
               values.Add(path.value(j));
            stats.add(values);
            EuropeanPathPricer pathPricer = new EuropeanPathPricer(Option.Type.Call, 110, riskFreeTS.link.discount(timeHorizon));
            double price = pathPricer.value(path);
            optionStat.add(price, 1);
         }
         List<double> discountedExpectedStockValue = new List<double>();
         List<double> expectedForward = stats.mean();
         for (int i = 0; i < steps; i++)
            discountedExpectedStockValue.Add(expectedForward[i] * riskFreeTS.link.discount(timeGrid[i]));

         for (int i = 0; i < steps; i++)
            Console.WriteLine("indice: {0}\t forward: {1}\t discounted: {2}\n", i, expectedForward[i], discountedExpectedStockValue[i]);
         TimeSpan span = DateTime.Now - now;
         Console.WriteLine("exec milliseconds: " + span.TotalMilliseconds);
         Console.ReadLine();
         
      }

      static void test()
      {
         #region Générateur d'uniformes
         RandomSequenceGenerator<MersenneTwisterUniformRng> URNG = new RandomSequenceGenerator<MersenneTwisterUniformRng>(250, 34);

         #endregion

         #region Générateur de loi normales par inversion de la fonction de répartition
         InverseCumulativeRsg<RandomSequenceGenerator<MersenneTwisterUniformRng>, InverseCumulativeNormal> GSG =
            new InverseCumulativeRsg<RandomSequenceGenerator<MersenneTwisterUniformRng>, InverseCumulativeNormal>(URNG, new InverseCumulativeNormal());
         for (uint i = 0; i < 4; i++)
         {
            Sample<List<double>> generatedSample = GSG.nextSequence();
            string str = "";
            foreach (double d in generatedSample.value)
               str += d + "\t";
            str += "\n";
            Console.WriteLine(str);
         }

         Console.ReadLine();
         Sample<List<double>> sample = GSG.nextSequence();
         #endregion

         #region Statistiques
         Statistics statistics = new Statistics();
         statistics.addSequence(sample.value, Enumerable.Repeat<double>(1, sample.value.Count).ToList());
         Console.WriteLine("test");
         double NormalSequenceVaR95 = statistics.valueAtRisk(0.975);
         Console.WriteLine(NormalSequenceVaR95);
         Console.ReadLine();
         #endregion
      }
   }
}
