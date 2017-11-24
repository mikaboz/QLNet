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
         #region Générateur d'uniformes
         RandomSequenceGenerator<MersenneTwisterUniformRng> URSG = new RandomSequenceGenerator<MersenneTwisterUniformRng>(250, 34);
         #endregion

         #region Générateur de loi normales par inversion de la fonction de répartition
         InverseCumulativeRsg<RandomSequenceGenerator<MersenneTwisterUniformRng>, InverseCumulativeNormal> NormalGenerator =
            new InverseCumulativeRsg<RandomSequenceGenerator<MersenneTwisterUniformRng>, InverseCumulativeNormal>(URSG, new InverseCumulativeNormal());
         for (uint i = 0; i < 4; i++)
         {
            Sample<List<double>> generatedSample = NormalGenerator.nextSequence();
            string str = "";
            foreach (double d in generatedSample.value)
               str += d + "\t";
            str += "\n";
            Console.WriteLine(str);
         }

         Console.ReadLine();
         Sample<List<double>> sample = NormalGenerator.nextSequence();
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
