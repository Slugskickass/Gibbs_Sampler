import org.jetbrains.kotlin.gradle.tasks.KotlinCompile

plugins {
    kotlin("jvm") version "1.4.0-rc"
}
group = "me.ashley"
version = "1.0-SNAPSHOT"

repositories {
    jcenter()
    mavenCentral()
    maven("https://dl.bintray.com/mipt-npm/scientifik")
    maven {
        url = uri("https://dl.bintray.com/kotlin/kotlin-eap")
    }
}
dependencies {
    testImplementation(kotlin("test-junit"))
    implementation("org.apache.commons:commons-csv:1.7")
    api("scientifik:kmath-core:${"0.1.3"}")
    implementation("org.apache.commons:commons-math3:3.6.1")
    implementation("no.tornado:tornadofx:1.7.20")
}
tasks.withType<KotlinCompile>() {
    kotlinOptions.jvmTarget = "1.8"
}