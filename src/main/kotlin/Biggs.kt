import org.apache.commons.math3.distribution.BinomialDistribution
import org.apache.commons.math3.distribution.BetaDistribution
import org.apache.commons.math3.distribution.NormalDistribution
import org.apache.commons.math3.distribution.GammaDistribution
import org.apache.commons.math3.stat.StatUtils
import org.apache.commons.math3.stat.descriptive.rank.Percentile
import kotlin.math.*
import tornadofx.*

fun nmean(Ds: DoubleArray, Dn: DoubleArray, cls: Double, cln: Double, ns: Int, nn: Int, pms: Pair<Double, Double>, pmn: Pair<Double, Double>): Pair<Double, Double> {
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

fun testmean(data: DoubleArray, prms: Pair<Double,Double>, prmn: Pair<Double,Double>): Pair<Double, Double> {
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
    var b1a = Dn.map{(it - cmn).pow(2)}
    var b1 = b1a.sum() / (2 + pln.second)
    //lan = np.random.gamma(shape=a1, scale=1 / b1)
    val lan = GammaDistribution(a1, 1/b1).sample()



    a1 = ns / 2 + pls.first
    //b1 = np.sum(np.square(Ds - cms * np.ones(ns))) / 2 + pls[1]
    val b1b = Ds.map{(it - cms).pow(2)}
    b1 = b1b.sum() / (2 + pls.second)
    val las = GammaDistribution(a1, 1 / b1).sample()

    return Pair(las, lan)
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
    var output = x.map { 0.5 * (ln(l) - ln(2 * 3.1415) - (it - m).pow(2)) }
    return output.toDoubleArray()
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

fun testarray(arraysee: DoubleArray, times: Int){
    for (i in 0 until times) {
        print(arraysee[i])
        print(' ')
    }
    println()
}

fun McMix(M: Int, X: DoubleArray, pms: Pair<Double, Double>, pmn: Pair<Double, Double>, pls: Pair<Double, Double>, pln: Pair<Double, Double>, pp: Pair<Double, Double>): Biggs {
    val N =M
    val data1 = Biggs("Miguel", N)
    var Ms = DoubleArray(N)
    var Mn = DoubleArray(N)
    var Ls = DoubleArray(N)
    var Ln = DoubleArray(N)
    var Pr = DoubleArray(N)
    var Id = DoubleArray(X.size)
    val percent = definePercentile(95.0)
    val threshold = percent.evaluate(X)

    val cuS = X.filter { it >= threshold }.toDoubleArray()
    val cuN = X.filter { it < threshold }.toDoubleArray()

    val Ns = cuS.size
    val Nn = cuN.size

    var cms = StatUtils.mean(cuS)
    var cmn = StatUtils.mean(cuN)

    var cls = 1 / StatUtils.variance(cuS)
    var cln = 1 / StatUtils.variance(cuN)

    val cPG = BetaDistribution(Ns + pp.first, Nn + pp.second)


    for (I in 0 until M) {

        val cp = cPG.sample()
        val CD = nmean(cuS, cuN, cls, cln, Ns, Nn, pms, pmn)
        cms = CD.first
        cmn = CD.second
        val CL = nprec(cuS, cuN, cms, cmn, Ns, Nn, pls, pln)
        cls = CL.first
        cln = CL.second

        val aux1 = gauspdf(X, cms, cls)
        val aux2 = gauspdf(X, cmn, cln)
        var aux = aux1.zip(aux2).map {it.first-it.second }.toDoubleArray()
        val pi = aux.map { 1 + ((1 - cp) / (cp * exp(it))) }.toDoubleArray()


        var cuSn = pi.map { BinomialDistribution(1, 1/it).sample().toInt() }
        val cuNn = cuSn.map{ abs(it-1)}

        val Ns = cuSn.sum()
        val Nn = cuNn.sum()

        // While loop
  //      Xs = X[np.where(cuS)]
  //      Xn = X[np.where(cuN)]
        data1.Ms[I] = cms
        data1.Mn[I] = cmn
        data1.Ls[I] = cls
        data1.Ln[I] = cln
        data1.Pr[I] = cp
        data1.Id += cuS
    }


return data1
}

data class Biggs(val name: String, val M: Int){
    var Ms = DoubleArray(M)
    var Mn = DoubleArray(M)
    var Ls = DoubleArray(M)
    var Ln = DoubleArray(M)
    var Pr = DoubleArray(M)
    var Id = DoubleArray(M)
}

class MyApp: App(MyView::class)

fun main() {

// Generate some test data
    val testdata  = generatetestdata(2000, 0.1, 1200, 0.005, 800, 0.01)
//
    val MM = 50                                           // Loops
//    val prms = Pair(1500, 10)               // Mean signal
//    val prmn = Pair(500, 1)                // mean noise
    val prms = Pair(definePercentile(75.0).evaluate(testdata), 0.1)               // Mean signal
    val prmn = Pair(definePercentile(25.0).evaluate(testdata), 0.1)                  // mean noise
    val prl = Pair(2.0, 0.6)        // Precisson
    val prp = Pair(0.5, 0.5)       // Precisson

//    val out = testmean(testdata, prms, prmn)
//    val ter =testnprec(testdata,out)

    val hold = McMix(MM,testdata, prms, prmn,prl,prl,prp)
    testarray(hold.Ms, MM)
//    launch<MyApp>()

}