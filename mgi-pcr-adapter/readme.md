### env install
```
sudo apt-get update
sudo apt-get install openjdk-8-jdk
sudo apt install maven
```
### change jdk version
```
sudo update-alternatives --config java

java --verison
```
### set maven java_home
```
ls /usr/lib/jvm/

vim ~/.mavenrc
export JAVA_HOME=/usr/lib/jvm/java-8-openjdk-amd64

mvn -v
```
###  modify the properties of spring/dockerfile
```
spring..url=...
SPRING_PCR_MYSQL_URL=...
```

