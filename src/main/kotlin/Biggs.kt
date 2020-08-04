import kotlin.math.ln
import kotlin.math.pow
import org.apache.commons.math3.distribution.BinomialDistribution
import org.apache.commons.math3.distribution.BetaDistribution
import org.apache.commons.math3.distribution.NormalDistribution
import org.apache.commons.math3.distribution.GammaDistribution
import org.apache.commons.math3.stat.StatUtils
import org.apache.commons.math3.stat.descriptive.rank.Percentile
import kotlin.math.sqrt

fun nmean(Ds: DoubleArray, Dn: DoubleArray, cls: Double, cln: Double, ns: Int, nn: Int, pms: Pair<Double, Double>, pmn: Pair<Double, Double>): Pair<Double, Double> {
    //p1 = pmn[1] + cln * nn
    var p1 = pmn.second + cln * nn

    //m1 = (pmn[0] * pmn[1] + cln * np.sum(Dn)) / p1
    var m1 = (pmn.first * pmn.second + cln * Dn.sum()) / p1

    //mn = m1 + np.random.standard_normal(1) / np.sqrt(p1)
    val mn = m1 + ( NormalDistribution().sample() / sqrt(p1))

    //p1 = pms[1] + cls * ns
    p1 = pmn.second + cls * ns

    //m1 = (pms[0] * pms[1] + cls * np.sum(Ds)) / p1
    m1 = (pms.first * pms.second + cls * Dn.sum()) / p1

    //ms = m1 + np.random.standard_normal(1) / np.sqrt(p1)
    val ms = m1 + ( NormalDistribution().sample() / sqrt(p1))

    //return ms, mn
    return Pair(ms,ms)
}

fun nprec(Ds: DoubleArray, Dn: DoubleArray, cms: Double, cmn: Double, ns: Int, nn: Int, pls: Pair<Double, Double>, pln: Pair<Double, Double>): Pair<Double, Double> {
    //a1 = nn / 2 + pln[0]
    val a1 = nn / 2 + pln.first

    //b1 = np.sum(np.square(Dn - cmn * np.ones(nn))) / 2 + pln[1]
    var counter = 0
    var b1a = DoubleArray(Dn.size)
    Dn.forEach {
        b1a[counter] = (it - cmn) * (it - cmn)
     counter  =+ 1
    }
    val b1 = b1a.sum() / (2 + pln.second)
    //print(a1, b1)

    //lan = np.random.gamma(shape=a1, scale=1 / b1)
    val lan = GammaDistribution(a1, 1/b1).sample()


    return Pair(lan, 0.0)
}
//def nprec(Ds, Dn, cms, cmn, ns, nn, pls, pln):



//lan = np.random.gamma(shape=a1, scale=1 / b1)
//# lan = np.random.gamma(shape=a1, scale=b1)
//
//a1 = ns / 2 + pls[0]
//b1 = np.sum(np.square(Ds - cms * np.ones(ns))) / 2 + pls[1]
//print(a1, b1)
//las = np.random.gamma(shape=a1, scale=1 / b1)
//# las = np.random.gamma(shape=a1, scale=b1)
//
//return las, lan



fun gauspdf(x: DoubleArray, m: DoubleArray, l: Double): DoubleArray {
    var counter: Int =0
    var output = DoubleArray(x.size)
    x.forEach {
        output[counter] = 0.5 * (ln(l) - ln(2 * 3.1415) - (it - m[counter]).pow(2))
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

fun McMix(M: Int, X: DoubleArray, pms: Pair<Int, Int>, pmn:Pair<Int, Int>, pls: Pair<Double, Double>, pln: Pair<Double, Double>, pp: Pair<Double, Double>): List<Double> {
    val N = X.size
    var Ms = DoubleArray(N)
    var Mn = DoubleArray(N)
    var Ls = DoubleArray(N)
    var Ln = DoubleArray(N)
    var Pr = DoubleArray(N)
    var Id = DoubleArray(N)
    val percent = definePercentile(75.0)
    val threshold = percent.evaluate(X)

    val cuS = X.filter { it >= threshold }
    val cuN = X.filter { it < threshold }

    val Ns = cuS.size
    val Nn = cuN.size

    val cls = 1 / StatUtils.variance(cuS.toDoubleArray())
    val cln = 1 / StatUtils.variance(cuN.toDoubleArray())

    val cP = BetaDistribution(Ns + pp.first, Nn + pp.second)

//    Xs       is cuS
//    Xn       is cuN
//    cls      value of 1 / variance of greater than terms
//    cln      value of 1 / variance of less than terms
//    Ns       number of greater than in array
//    Nn       number of less than in array
//    pms      The pair with
//    pmn      The pair
//     nmean(Ds, Dn, cls, cln, ns, nn, pms, pmn)


    for (i in 0 until M) {

        println(cP.sample())

    }

return cuS
}

fun main() {
    val MM = 5                              // Loops
    val nn = 15000
    val prms = Pair(1500, 10)
    val prmn = Pair(500, 1)
    val prl = Pair(2.0, 0.6)   // Means
    val prp = Pair(0.5, 0.5)    // Prec

    val testdata  = generatetestdata(15000, 0.5, 1200, 0.005, 800, 0.001)

    val test = McMix(MM, testdata,prms, prmn, prl, prl, prp)
    println(test.size)

}