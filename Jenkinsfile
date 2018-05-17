pipeline {
    agent none
    environment {
        PETSC_DIR = '/usr/lib/petscdir/3.8.3'
    }
    stages {
        stage('Building') {
            agent { label 'azure-linux' }
            steps {
                script {
                    def version = sh script: 'git rev-parse --short HEAD | tr -d "\n"', returnStdout: true
                    def imageName = "icokpy"
                    def imageTag = "${env.BUILD_NUMBER}_${version}"
                    def fluidityImage = docker.build("${imageName}:${imageTag}")

                    fluidityImage.inside() {
                        sh 'make unittest'
                    }
                }
            }
        }
    }
}
