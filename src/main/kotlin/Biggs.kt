import javafx.scene.chart.Chart
import javafx.scene.layout.VBox
import kotlin.math.ln
import kotlin.math.pow
import org.apache.commons.math3.distribution.BinomialDistribution
import org.apache.commons.math3.distribution.BetaDistribution
import org.apache.commons.math3.distribution.NormalDistribution
import org.apache.commons.math3.distribution.GammaDistribution
import org.apache.commons.math3.stat.StatUtils
import org.apache.commons.math3.stat.descriptive.rank.Percentile
import kotlin.math.sqrt
import kotlin.math.exp

fun nmean(Ds: DoubleArray, Dn: DoubleArray, cls: Double, cln: Double, ns: Int, nn: Int, pms: Pair<Int, Int>, pmn: Pair<Int, Int>): Pair<Double, Double> {
    //p1 = pmn[1] + cln * nn
    var p1 = pmn.second + cln * nn
    //m1 = (pmn[0] * pmn[1] + cln * np.sum(Dn)) / p1
    var m1 = (pmn.first * pmn.second + cln * Dn.sum()) / p1
    //mn = m1 + np.random.standard_normal(1) / np.sqrt(p1)
    val mn = m1 + ( NormalDistribution().sample() / sqrt(p1))



    //p1 = pms[1] + cls * ns
    p1 = pms.second + cls * ns
    //m1 = (pms[0] * pms[1] + cls * np.sum(Ds)) / p1
    m1 = (pms.first * pms.second + cls * Ds.sum()) / p1
    //ms = m1 + np.random.standard_normal(1) / np.sqrt(p1)
    val ms = m1 + ( NormalDistribution().sample() / sqrt(p1))

    return Pair(ms,mn)
}

fun testmean(data: DoubleArray, prms: Pair<Int,Int>, prmn: Pair<Int,Int>): Pair<Double, Double> {
    val percent = definePercentile(75.0)
    val threshold = percent.evaluate(data)

    val cuS = data.filter { it >= threshold }.toDoubleArray()
    val cuN = data.filter { it < threshold }.toDoubleArray()

    val cls = 1 / StatUtils.variance(cuS)
    val cln = 1 / StatUtils.variance(cuN)

    val Ns = cuS.size
    val Nn = cuN.size

    val outdata = nmean(cuS, cuN, cls, cln, Ns, Nn, prms, prmn)
    //println(outdata.first)
    //println(outdata.second)

    return outdata
}

fun nprec(Ds: DoubleArray, Dn: DoubleArray, cms: Double, cmn: Double, ns: Int, nn: Int, pls: Pair<Double, Double>, pln: Pair<Double, Double>): Pair<Double, Double> {
    //a1 = nn / 2 + pln[0]
    var a1 = nn / 2 + pln.first
    //b1 = np.sum(np.square(Dn - cmn * np.ones(nn))) / 2 + pln[1]
    var counter = 0
    var b1a = DoubleArray(Dn.size)
    Dn.forEach {
        b1a[counter] = (it - cmn).pow(2)
        counter += 1
    }
    var b1 = b1a.sum() / (2 + pln.second)
    //lan = np.random.gamma(shape=a1, scale=1 / b1)
    val lan = GammaDistribution(a1, 1/b1).sample()





    a1 = ns / 2 + pls.first

    //b1 = np.sum(np.square(Ds - cms * np.ones(ns))) / 2 + pls[1]
    counter = 0
    b1a = DoubleArray(Dn.size)
    Ds.forEach {
        b1a[counter] = (it - cms).pow(2)
        counter  += 1
    }


    b1 = b1a.sum() / (2 + pls.second)

    val las = GammaDistribution(a1, 1 / b1).sample()

    return Pair(lan, las)
}

fun testnprec(data: DoubleArray, out: Pair<Double,Double>, ): Pair<Double, Double> {

    val percent = definePercentile(75.0)
    val threshold = percent.evaluate(data)
    val cuS = data.filter { it >= threshold }.toDoubleArray()
    val cuN = data.filter { it < threshold }.toDoubleArray()

    val cls = 1 / StatUtils.variance(cuS)
    val cln = 1 / StatUtils.variance(cuN)

    val Ns = cuS.size
    val Nn = cuN.size

    val prp = Pair(0.5, 0.5)

    val ter =nprec(cuS, cuN, out.first, out.second, Ns, Nn,prp, prp)
    return ter
}

fun gauspdf(x: DoubleArray, m: Double, l: Double): DoubleArray {
    var counter: Int =0
    var output = DoubleArray(x.size)
    x.forEach {
        output[counter] = 0.5 * (ln(l) - ln(2 * 3.1415) - (it - m).pow(2))
        counter += 1
    }
    return output
}

fun playBinomialDistribution(trails: Int, p: Double, runs: Int): DoubleArray {
    val outer = BinomialDistribution(trails, p)
    var tempting = DoubleArray(runs)
    for (i in 0 until runs) {
        tempting[i] = outer.sample().toDouble()
    }
    return tempting
}

fun definePercentile(P: Double): Percentile {
    val out =(Percentile(P))
    return out
}

fun generatetestdata(points: Int, pp: Double, mus: Int, las: Double, mun: Int, lan: Double): DoubleArray{
    val newt = playBinomialDistribution(1,pp, points)
    val newg =NormalDistribution()
    var mall = DoubleArray(points)
    var laall = DoubleArray(points)
    var counter = 0
    newt.forEach {
        mall[counter] = (mus * it) + (mun * (1-it))
        laall[counter] = (las * it) + (lan * (1-it))
        counter += 1}

    var XX = DoubleArray(points)
    counter = 0
    mall.forEach {
        val tempdata = newg.sample()
        XX[counter] = (it + (tempdata/ sqrt(laall[counter])))
        counter += 1
    }
    return XX
}

fun McMix(M: Int, X: DoubleArray, pms: Pair<Int, Int>, pmn:Pair<Int, Int>, pls: Pair<Double, Double>, pln: Pair<Double, Double>, pp: Pair<Double, Double>): Double {
    val N = X.size
    var Ms = DoubleArray(N)
    var Mn = DoubleArray(N)
    var Ls = DoubleArray(N)
    var Ln = DoubleArray(N)
    var Pr = DoubleArray(N)
    var Id = DoubleArray(N)
    val percent = definePercentile(75.0)
    val threshold = percent.evaluate(X)

    val cuS = X.filter { it >= threshold }.toDoubleArray()
    val cuN = X.filter { it < threshold }.toDoubleArray()

    val Ns = cuS.size
    val Nn = cuN.size

    val cms = StatUtils.mean(cuS)
    val cmn = StatUtils.mean(cuN)

    val cls = 1 / StatUtils.variance(cuS)
    val cln = 1 / StatUtils.variance(cuN)

    val cPG = BetaDistribution(Ns + pp.first, Nn + pp.second)


    for (i in 0 until M) {

        val cp = cPG.sample()
        val CD = nmean(cuS, cuN, cls, cln, Ns, Nn, pms, pmn)
        val CL = nprec(cuS, cuN, cms, cmn, Ns, Nn, pls, pln)

        val aux1 = gauspdf(X, cms, cls)
        val aux2 = gauspdf(X, cmn, cln)
        var aux = aux1.zip(aux2).map {it.first-it.second }.toDoubleArray()


 //       val pi = 1 + ((1 - cp) / cp * exp(aux))

        var pie =DoubleArray(aux.size)
        var counter = 0
        aux.forEach {pie[counter] = 1 + ((1 - cp) / cp * exp(it))
            counter =+1
        }


        val pi = pie.filter { it.isFinite() }


       var cuSn = IntArray(pi.size)
       counter = 0
       pi.forEach {
//           cuSn[counter] = BinomialDistribution(1, 1 / it).sample()
           cuSn[counter] = BinomialDistribution(1, it).sample()
       counter =+1
       }

        cuSn.forEach { println(it) }

//        cuNn = cuSn.forEach { not(it) }
//        Ns = cuSn.sum()
//        Nn = cuNn.sum()

    }

return 0.0
}



fun main() {
// Generate some test data
    val testdata  = generatetestdata(2000, 0.1, 1200, 0.005, 800, 0.01)
//
    val MM = 5                                         // Loops
    val nn = 2000                                      // Points
    val prms = Pair(1500, 10)              // Mean signal
    val prmn = Pair(500, 1)                // mean noise
    val prl = Pair(2.0, 0.6)        // Precisson
    val prp = Pair(0.5, 0.5)        // Precisson

    val out = testmean(testdata, prms, prmn)
    val ter =testnprec(testdata,out)

}