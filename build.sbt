import org.scalajs.sbtplugin.ScalaJSPlugin.autoImport.*
import sbtcrossproject.CrossPlugin.autoImport.*
import scalajscrossproject.ScalaJSCrossPlugin.autoImport.*

ThisBuild / organization := "io.github.canardlapin"
ThisBuild / scalaVersion := "3.4.2"
ThisBuild / version      := "0.1.0-SNAPSHOT"
ThisBuild / versionScheme := Some("early-semver")
ThisBuild / homepage := Some(url("https://github.com/canardlapin/multivar"))
ThisBuild / licenses := Seq(
  "Apache-2.0" -> url("https://www.apache.org/licenses/LICENSE-2.0.txt")
)
ThisBuild / scmInfo := Some(
  ScmInfo(
    url("https://github.com/canardlapin/multivar"),
    "scm:git:https://github.com/canardlapin/multivar.git",
    Some("scm:git:git@github.com:canardlapin/multivar.git")
  )
)
ThisBuild / developers := List(
  Developer(
    "canardlapin",
    "canardlapin",
    "307091466+canardlapin@users.noreply.github.com",
    url("https://github.com/canardlapin")
  )
)
ThisBuild / pomIncludeRepository := (_ => false)
ThisBuild / publishMavenStyle := true

// Immutable source dependency containing the shared matrix, operator, solver,
// spectral, and first-order optimization substrate.
lazy val galeRevision = "d55fe2f97196a76ab7879e1a12f1e92403aeba06"
lazy val galeBuild =
  uri(s"https://github.com/canardlapin/gale.git#$galeRevision")
lazy val galeCoreJVM = ProjectRef(galeBuild, "coreJVM")
lazy val galeCoreJS  = ProjectRef(galeBuild, "coreJS")

// Multivar retains ScalaFIM's deliberate 3.4.2 consumer line. Gale's portable
// sources are cross-checked on that compiler even though Gale itself publishes
// its current baseline from Scala 3.7.4.
galeCoreJVM / scalaVersion := "3.4.2"
galeCoreJS / scalaVersion  := "3.4.2"

lazy val commonSettings = Seq(
  scalacOptions ++= Seq(
    "-deprecation",
    "-feature",
    "-unchecked",
    "-Xmax-inlines:64"
  ),
  Test / fork := false,
  libraryDependencies += "org.scalameta" %%% "munit" % "1.2.1" % Test
)

lazy val jsSettings = Seq(
  scalaJSLinkerConfig ~= (_.withModuleKind(ModuleKind.CommonJSModule)),
  Test / jsEnv := new org.scalajs.jsenv.nodejs.NodeJSEnv()
)

lazy val core =
  crossProject(JSPlatform, JVMPlatform)
    .crossType(CrossType.Full)
    .in(file("modules/core"))
    .settings(commonSettings)
    .settings(
      name := "multivar-core",
      description := "Typed, evidence-bearing multivariate analysis for Scala."
    )
    .jvmConfigure(_.dependsOn(galeCoreJVM))
    .jsConfigure(_.dependsOn(galeCoreJS))
    .jsSettings(jsSettings)

lazy val coreJS  = core.js
lazy val coreJVM = core.jvm

lazy val ir =
  crossProject(JSPlatform, JVMPlatform)
    .crossType(CrossType.Full)
    .in(file("modules/ir"))
    .dependsOn(core)
    .settings(commonSettings)
    .settings(
      name := "multivar-ir",
      description := "Language-neutral schemas and codecs for multivar programs and evidence."
    )
    .jvmConfigure(_.dependsOn(galeCoreJVM))
    .jsConfigure(_.dependsOn(galeCoreJS))
    .jsSettings(jsSettings)

lazy val irJS  = ir.js
lazy val irJVM = ir.jvm

lazy val root =
  project
    .in(file("."))
    .aggregate(coreJVM, coreJS, irJVM, irJS)
    .settings(
      name := "multivar",
      publish / skip := true
    )

addCommandAlias("compileAll", ";coreJVM/compile;coreJS/compile;irJVM/compile;irJS/compile")
addCommandAlias("testAll", ";coreJVM/test;coreJS/test;irJVM/test;irJS/test")
