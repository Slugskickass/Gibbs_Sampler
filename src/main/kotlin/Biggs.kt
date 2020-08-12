import org.apache.commons.math3.distribution.BinomialDistribution
import org.apache.commons.math3.distribution.BetaDistribution
import org.apache.commons.math3.distribution.NormalDistribution
import org.apache.commons.math3.distribution.GammaDistribution
import org.apache.commons.math3.stat.StatUtils
import org.apache.commons.math3.stat.descriptive.rank.Percentile
import kotlin.math.*
import tornadofx.*

// TODO: Make Compile
// TODO: newMeans data class
// TODO: newMeans extract function
// TODO: more data classes (variance / precision)

fun newMeans(signalSamples: DoubleArray, noiseSamples: DoubleArray, signalPrior: Prior, noisePrior: Prior): Pair<Double, Double> {
    val signalPrecision = 1.0 / StatUtils.variance(signalSamples)
    val noisePrecision = 1.0 / StatUtils.variance(noiseSamples)

    val updatedNoisePrecision =
        noisePrior.precision + noisePrecision * noiseSamples.size
    val updatedNoiseMean =
        (noisePrior.mean * noisePrior.precision + noisePrecision * noiseSamples.sum()) / updatedNoisePrecision
    val noiseMean =
        NormalDistribution(updatedNoiseMean, 1.0 / sqrt(updatedNoisePrecision)).sample()

    val updatedSignalPrecision =
        signalPrior.precision + signalPrecision * signalSamples.size
    val updatedSignalMean =
        (signalPrior.mean * signalPrior.precision + signalPrecision * signalSamples.sum()) / updatedSignalPrecision
    val sampleMean =
        NormalDistribution(updatedSignalMean, 1.0 / sqrt(updatedSignalPrecision)).sample()

    return Pair(max(noiseMean, sampleMean), min(noiseMean, sampleMean))
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

data class Probability(val p : Double) {
    init { assert(p in 0.0..1.0) }

    val np get() = 1.0 - p
}

data class Prior(val mean : Double, val precision : Double)

data class McMixConfig(
    val numberOfIterations : Int,
    val thresholdPercentile : Percentile,
    val initialGuess : Probability,
    val signalPrior : Prior,
    val noisePrior : Prior
) {
    constructor(
        numberOfIterations : Int,
        thresholdPercentile : Percentile,
        initialGuess : Probability,
        data : DoubleArray,
        signalPriorPercentile : Percentile = definePercentile(75.0),
        signalPriorPrecision : Double = 0.1,
        noisePriorPercentile : Percentile = definePercentile(25.0),
        noisePriorPrecision : Double = 0.1,
    ) : this (
        numberOfIterations,
        thresholdPercentile,
        initialGuess,
        Prior(signalPriorPercentile.evaluate(data), signalPriorPrecision),
        Prior(noisePriorPercentile.evaluate(data), noisePriorPrecision)
    )
}

fun runMcMix(config : McMixConfig, input : DoubleArray) : Biggs {
    val output = Biggs("Miguel", config.numberOfIterations)
    val threshold = config.thresholdPercentile.evaluate(input)

    val geqThreshold = input.filter { it >= threshold }.toDoubleArray()
    val ltThreshold = input.filter { it < threshold }.toDoubleArray()

    val betaDistribution = BetaDistribution(
        geqThreshold.size + config.initialGuess.p,
        ltThreshold.size + config.initialGuess.np
    )

    for (I in 0 until config.numberOfIterations) {
        val betaSample = betaDistribution.sample()
        val CD = newMeans(geqThreshold, ltThreshold, config.signalPrior, config.noisePrior)
        val cms = CD.first
        val cmn = CD.second
        val CL = nprec(cuS, cuN, cms, cmn, Ns, Nn, pls, pln)
        cls = CL.first
        cln = CL.second

        val aux1 = gauspdf(X, cms, cls)
        val aux2 = gauspdf(X, cmn, cln)
        var aux = aux1.zip(aux2).map {it.first-it.second }.toDoubleArray()
        val pi = aux.map { 1 + ((1 - betaSample) / (betaSample * exp(it))) }.toDoubleArray()


        var cuS = pi.map { BinomialDistribution(1, 1/it).sample().toInt() }
        val cuN = cuS.map{ abs(it-1)}

        val Ns = cuS.sum()
        val Nn = cuN.sum()

        // While loop
        //      Xs = X[np.where(cuS)]
        //      Xn = X[np.where(cuN)]
        output.meanSignal[I] = cms
        output.meanNoise[I] = cmn
        output.signalVariance[I] = cls
        output.noiseVariance[I] = cln
        output.precision[I] = betaSample
        //      data1.Id += cuS
        for (i in 0 until output.pSignal.size) {
            output.pSignal[i] = cuS[i].toDouble() + output.pSignal[i]
        }
    }


    return output
}

fun McMix(M: Int, X: DoubleArray, pms: Pair<Double, Double>, pmn: Pair<Double, Double>, pls: Pair<Double, Double>, pln: Pair<Double, Double>, pp: Pair<Double, Double>): Biggs {
    val N = M
    var Ms = DoubleArray(N)
    var Mn = DoubleArray(N)
    var Ls = DoubleArray(N)
    var Ln = DoubleArray(N)
    var Pr = DoubleArray(N)
    var Id = DoubleArray(X.size)
    val threshold = definePercentile(95.0).evaluate(X)

    val cuS = X.filter { it >= threshold }.toDoubleArray()
    val cuN = X.filter { it < threshold }.toDoubleArray()

    val Ns = cuS.size
    val Nn = cuN.size

//    var cms = StatUtils.mean(cuS)
//    var cmn = StatUtils.mean(cuN)

    var cls = 1 / StatUtils.variance(cuS)
    var cln = 1 / StatUtils.variance(cuN)

    val cPG = BetaDistribution(Ns + pp.first, Nn + pp.second)
}

data class Biggs(val name: String, val M: Int) {
    var meanSignal = DoubleArray(M)
    var meanNoise = DoubleArray(M)
    var signalVariance = DoubleArray(M)
    var noiseVariance = DoubleArray(M)
    var precision = DoubleArray(M)
    var pSignal = DoubleArray(M)
}

class MyApp: App(MyView::class)

fun main() {

// Generate some test data
    val testdata  = generatetestdata(2000, 0.3, 1500, 0.1, 500, 0.05)

    val MM = 200                                           // Loops
    val prms = Pair(definePercentile(75.0).evaluate(testdata), 0.1)               // Mean signal
    val prmn = Pair(definePercentile(25.0).evaluate(testdata), 0.1)                  // mean noise
    val prl = Pair(2.0, 0.6)        // Precisson
    val prp = Pair(0.5, 0.5)       // Precisson

//    val out = testmean(testdata, prms, prmn)
//    val ter =testnprec(testdata,out)

    val hold = McMix(MM,testdata, prms, prmn,prl,prl,prp)

    println(hold.meanSignal[MM-1])
    println(hold.meanNoise[MM-1])
    testarray(hold.pSignal, 30)
}