import javafx.scene.chart.NumberAxis
import tornadofx.*
import java.awt.Button


class MyView : View("Gibbs Sampler"){
    override val root = borderpane() {
        top<TopView>()
        bottom<BottomView>()
    }
    class TopView: View() {
            override val root = hbox{
            button {
                this.text = "Generate Data"
                action { val testdata  = generatetestdata(2000, 0.1, 1200, 0.005, 800, 0.01) }
            }
                button {
                    this.text = "Show Data"

            }



                button {
                    this.text = "Run Gibbs"
                    action { println("hello") }
                }
                button {
                    this.text = "Show Gibbs"
                    action { println("hello") }
                }
        }
        }

    class BottomView: View() {
        override val root = hbox{
            scatterchart("Noise over time", NumberAxis(), NumberAxis()) {
                series("Noise") {
                    var count = 0
                    //               testdata.forEach{
                    //                   count +=1
                    //                   data(count, it)
                }
            }
        }

}}
class changedat{
    val testdata  = generatetestdata(2000, 0.1, 1200, 0.005, 800, 0.01)

}


//class MyView : View("My View") {
//    override val root = hbox {
//        val testdata  = generatetestdata(2000, 0.1, 1200, 0.005, 800, 0.01)
//        val MM = 50                                           // Loops
//        val nn = 2000                                       // Points
//        val prms = Pair(definePercentile(75.0).evaluate(testdata), 0.1)               // Mean signal
//        val prmn = Pair(definePercentile(25.0).evaluate(testdata), 0.1)                  // mean noise
//        val prl = Pair(2.0, 0.6)        // Precisson
//        val prp = Pair(0.5, 0.5)       // Precisson
//        val hold = McMix(MM,testdata, prms, prmn,prl,prl,prp)
//        val x_data = IntArray(MM, { i -> i })
//
//        scatterchart("Noise over time", NumberAxis(), NumberAxis()) {
//            series("Noise") {
//                var count = 0
//                hold.Mn.forEach{
//                    count +=1
//                    data(count, it)
//                }
//            }
//        }
//        scatterchart("Signal over time", NumberAxis(), NumberAxis()) {
//            series("Noise") {
//                var count = 0
//                hold.Ms.forEach {
//                    count += 1
//                    data(count, it)
//                }
//            }
//        }
//        button {
//        this.text = "test"
//        action { println("hello") }}
//    }
//}

