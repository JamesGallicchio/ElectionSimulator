import Stats.MultiNormal

import LinearAlgebra._

object ElectionSimulator {

  type Region = String
  type Candidate = String
  case class Poll(name: String, tag: Region, weight: Double,
                  results: Map[Candidate, Double])

  case class State(name: String, atLarge: Int, districts: Set[District])
  case class District(name: String, delegates: Int, pollTags: Map[Region, Double])

  def associatePolls(districts: Set[District], polls: Set[Poll]
                    ): Map[District, Map[Poll, Double]] =
    districts.map { d =>
      (d, d.pollTags.flatMap { case (reg, weight) =>
            polls.filter(_.tag == reg).map(p => (p, p.weight*weight)) })
    }.toMap

  def calculateDistrib(polls: Map[Poll, Double], candidates: Seq[Candidate])
                      (n: Int = candidates.size): MultiNormal[n.type] =
    MultiNormal.fromSamplesWeighted(n,
      polls.map { case (p: Poll, w) => (
          Vector(n, candidates.map(p.results.getOrElse(_, 0))).stochastic,
        w) }.toSeq
    )

  case class ElectionResults() {
    lazy val distDelegates: Map[District, Map[Candidate, Int]] = ???
    lazy val totalDelegates: Map[Candidate, Int] = ???
  }

  def simulateElection[N](districts: Map[District, MultiNormal[N]], candidates: Seq[Candidate]): ElectionResults =
    ElectionResults(districts.view.mapValues(v => (candidates zip v.sample.elems).toMap ).toMap)

  def main(args: Array[String]): Unit = {

  }
}
