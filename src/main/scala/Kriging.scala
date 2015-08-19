import scala.annotation.tailrec
import scala.collection.mutable.HashMap
import scala.collection.mutable.ArrayBuffer
import scala.io.Source

/**
 * Created by lenovo on 8/14/2015.
 */
object Kriging {
  def main(args: Array[String]) {
    val filePath = "data/ZoneA.dat"

    //create a matrix
    val rows: ArrayBuffer[Array[Double]] = new ArrayBuffer[Array[Double]]()
    for (line <- Source.fromFile(filePath).getLines().toList.drop(10)){
      rows.append(line.split(" ").slice(0,4).map(x => x.toDouble))
    }

    val bandwidth = 500
    val distMat: Array[Array[Double]] = distMatrix(rows.map(x => new Tuple2(x(0), x(1))).toArray)

    val lags = 0.until(10500, bandwidth).toList

    // TODO refactor into functional form
    // calculate semivariance on each lag
    // result: semivariogram
    val laggedSemivariance: List[Double] = lags.map { lag =>
      val semivariances: ArrayBuffer[Double] = new ArrayBuffer[Double]()
      for (i <- 0 until distMat.length) {
        for (j <- i+1 until distMat.length) {
          if (distMat(i)(j) >= lag-bandwidth && distMat(i)(j) <= lag+bandwidth) {
            val vars = Math.pow(rows(i)(3)-rows(j)(3), 2)
            semivariances.append(vars)
          }
        }
      }
      semivariances.sum / (2 * semivariances.length)
    }

    val semivariogram = Tuple2(lags, laggedSemivariance)

    val dataVariance = variance(rows.map{x => x(3)}.toList)
    getCovFunction(rows.map{x => x(3)}.toList, semivariogram, bandwidth)
  }

  // calculate distance matrix on spatial data
  def distMatrix(vals: Array[Tuple2[Double, Double]]): Array[Array[Double]] = {
    val distRow = new Array[Array[Double]](vals.length)
    for (i <- 0 until vals.length) {
      val distCol = new Array[Double](vals.length)
      for (j <- 0 until vals.length) {
        distCol(j) = Math.sqrt(Math.pow(vals(i)._1-vals(j)._1, 2)+Math.pow(vals(i)._2-vals(j)._2, 2))
      }
      distRow(i) = distCol
    }
    return distRow
  }

  def variance (vals: List[Double]): Double = {
    val mean = vals.sum / vals.length
    @tailrec
    def varTail(accu: Double, current: List[Double]): Double = {
      if (current.isEmpty) {
        return accu
      } else {
        return varTail(accu+Math.pow(current.head-mean, 2), current.tail)
      }
    }
    return (1.0 / vals.length) * varTail(0, vals)
  }

  def covariance (vals: List[Double], lag: Int, semivariogram: Tuple2[List[Int], List[Double]], bw: Int): Double = {
    val cov0 = variance(vals)
    return if (lag == 0) cov0 else cov0 - semivariogram._2(lag / bw)
  }

  // TODO spherical is the default model for the time being
  def getCovFunction (vals: List[Double], semivariogram: Tuple2[List[Int], List[Double]], bandwidth: Int) = {
    // calculate the sill
    val cov0 = covariance(vals, semivariogram._1(0), semivariogram, bandwidth)
    // calculate optimal parameters (practical range)
    val param = {
      val start = semivariogram._1(1)
      val end = semivariogram._1.last
      val meshSize = 1000
      val spacedList: List[Double] = {
        val space = (end-start).toDouble / (meshSize-1).toDouble
        List.tabulate(meshSize){n => start+(n*space)}
      }

      val mse: Array[Double] = Array.fill(meshSize){0}
      for (i <- 0 until spacedList.length) {
        mse(i) = {
          val errors = spherical(semivariogram._1, spacedList(i), cov0)
            .zip(semivariogram._2)
            .map { yResultTuple =>
              val y = yResultTuple._2
              val result = yResultTuple._1
              Math.pow( y - result, 2 )
            }
          errors.sum / errors.length
        }
      }
      spacedList( mse.zipWithIndex.min._2 )
    }
    println(param)
  }

  def spherical (lags: List[Int], examinedRange: Double, sill: Double): List[Double] = {
    def computeSpherical (lag: Int): Double = {
      if (lag <= examinedRange) {
        return sill * ( 1.5 * (lag/examinedRange) - 0.5 * Math.pow( (lag/examinedRange), 3 ) )
      } else {
        return sill
      }
    }
    // compute spherical function for each lags
    return lags.map{lag =>
      computeSpherical(lag)
    }
  }
}