import org.scalajs.sbtplugin.ScalaJSPlugin.autoImport.*
import org.typelevel.sbt.TypelevelSitePlugin
import org.typelevel.sbt.TypelevelSitePlugin.autoImport.*
import sbtcrossproject.CrossPlugin.autoImport.*
import scalajscrossproject.ScalaJSCrossPlugin.autoImport.*
import laika.ast.Path.Root
import laika.helium.config.{HeliumIcon, IconLink}

ThisBuild / organization := "io.github.canardlapin"
ThisBuild / scalaVersion := "3.7.4"
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
ThisBuild / tlSitePublishBranch := None
ThisBuild / tlSitePublishTags := false

// Gale is a Maven dependency so published multivar POMs declare a resolvable
// coordinate rather than a Git ProjectRef. Until Gale is on Maven Central, CI
// and local builds install the pinned revision with tools/publish-gale-local.sh.
lazy val galeRevision = "d55fe2f97196a76ab7879e1a12f1e92403aeba06"
lazy val galeVersion  = s"1.0.0-${galeRevision.take(12)}"

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

lazy val mimaSettings = Seq(
  // No previous release yet: MiMa stays idle until 0.1.0 ships, then point
  // mimaPreviousArtifacts at that coordinate.
  mimaPreviousArtifacts := Set.empty,
  mimaFailOnNoPrevious := false
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
      description := "Typed, evidence-bearing multivariate analysis for Scala.",
      libraryDependencies += "io.github.canardlapin" %%% "gale-core" % galeVersion
    )
    .jvmSettings(mimaSettings)
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
    .jvmSettings(mimaSettings)
    .jsSettings(jsSettings)

lazy val irJS  = ir.js
lazy val irJVM = ir.jvm

lazy val stageDocsApi = taskKey[Unit](
  "Stage the generated JVM Scaladoc for the public modules inside the guide site."
)

lazy val docs =
  project
    .in(file("site"))
    .dependsOn(coreJVM, irJVM)
    .enablePlugins(TypelevelSitePlugin)
    .settings(
      name := "multivar-docs",
      publish / skip := true,
      // mdoc examples import gale.linalg.*; TypelevelSitePlugin's mdoc classpath
      // does not always surface transitive libraryDependencies of dependsOn
      // projects after the Phase 4 Maven pin, so pin Gale here explicitly.
      libraryDependencies += "io.github.canardlapin" %% "gale-core" % galeVersion,
      mdocIn := (ThisBuild / baseDirectory).value / "site-docs",
      laikaConfig := LaikaConfig.defaults.withRawContent,
      tlSiteHelium := tlSiteHelium.value
        .site
        .topNavigationBar(
          homeLink = IconLink.internal(Root / "index.md", HeliumIcon.home)
        )
        .site
        .mainNavigation(
          depth = 3,
          includePageSections = false
        ),
      stageDocsApi := {
        val destination = mdocOut.value / "api"
        IO.delete(destination)
        IO.copyDirectory(
          (coreJVM / Compile / doc).value,
          destination / "core"
        )
        IO.copyDirectory(
          (irJVM / Compile / doc).value,
          destination / "ir"
        )
      },
      tlSite := Def
        .sequential(
          mdoc.toTask(""),
          stageDocsApi,
          laikaSite
        )
        .value
    )

/** Consumer that depends only on publishedLocal Maven/Ivy artifacts.
  *
  * It must not `dependsOn` the in-repo modules. A green `smoke/compile` is
  * evidence that the published dependency graph resolves without the source
  * tree that produced it.
  */
lazy val smoke =
  project
    .in(file("modules/smoke"))
    .settings(
      name := "multivar-smoke",
      publish / skip := true,
      scalaVersion := (ThisBuild / scalaVersion).value,
      scalacOptions ++= Seq("-deprecation", "-feature", "-unchecked"),
      libraryDependencies ++= Seq(
        "io.github.canardlapin" %% "multivar-core" % version.value,
        "io.github.canardlapin" %% "multivar-ir" % version.value
      )
    )

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
addCommandAlias("docsCheck", "docs/tlSite")
addCommandAlias(
  "smokeCheck",
  ";coreJVM/publishLocal;irJVM/publishLocal;smoke/clean;smoke/compile"
)
addCommandAlias("mimaCheck", ";coreJVM/mimaReportBinaryIssues;irJVM/mimaReportBinaryIssues")
