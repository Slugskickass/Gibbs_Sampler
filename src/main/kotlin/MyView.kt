import javafx.scene.chart.NumberAxis
import tornadofx.*

class MyView : View("My View") {
    override val root = hbox {
        val testdata  = generatetestdata(2000, 0.1, 1200, 0.005, 800, 0.01)
        val MM = 50                                           // Loops
        val nn = 2000                                       // Points
        val prms = Pair(definePercentile(75.0).evaluate(testdata), 0.1)               // Mean signal
        val prmn = Pair(definePercentile(25.0).evaluate(testdata), 0.1)                  // mean noise
        val prl = Pair(2.0, 0.6)        // Precisson
        val prp = Pair(0.5, 0.5)       // Precisson
        val hold = McMix(MM,testdata, prms, prmn,prl,prl,prp)
        val x_data = IntArray(MM, { i -> i })

        scatterchart("Noise over time", NumberAxis(), NumberAxis()) {
            series("Noise") {
                var count = 0
                hold.Mn.forEach{
                    count +=1
                    data(count, it)
                }
            }
        }
        scatterchart("Signal over time", NumberAxis(), NumberAxis()) {
            series("Noise") {
                var count = 0
                hold.Ms.forEach {
                    count += 1
                    data(count, it)
                }
            }
        }
        button {
        this.text = "test"
        action { println("hello") }}
    }
}

