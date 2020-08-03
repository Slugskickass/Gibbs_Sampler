import kotlin.math.ln
import kotlin.math.pow
import org.apache.commons.math3.distribution.BinomialDistribution
import org.apache.commons.math3.distribution.NormalDistribution
import org.apache.commons.math3.stat.descriptive.rank.Percentile
import kotlin.math.sqrt

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
        println(XX[counter])
    }
    return XX
}

fun McMix(M: Int, X: DoubleArray, pms: Pair<Int, Int>, pmn:Pair<Int, Int>, pls: Pair<Double, Double>, pln: Pair<Double, Double>, pp: Pair<Double, Double>){
    val N = X.size
    var Ms = DoubleArray(N)
    var Mn = DoubleArray(N)
    var Ls = DoubleArray(N)
    var Ln = DoubleArray(N)
    var Pr = DoubleArray(N)
    var Id = DoubleArray(N)
    val Threshold = definePercentile(0.74)

return
}

fun main() {
    val MM = 5                              // Loops
    val nn = 15000
    val prms = Pair(1500, 10)
    val prmn = Pair(500, 1)
    val prl = Pair(2, 0.6)   // Means
    val prp = Pair(1 / 2, 1 / 2)    // Prec

    val testdata  = generatetestdata(15000, 0.5, 1200, 0.005, 800, 0.001)
    print(testdata[100])
}