pipeline {
    agent none
    environment {
        PETSC_DIR = '/usr/lib/petscdir/3.8.3'
    }
    stages {
        stage('Building') {
          agent { label 'azure-linux' }
            steps {
                sh 'true'
            }
        }
    }
}
