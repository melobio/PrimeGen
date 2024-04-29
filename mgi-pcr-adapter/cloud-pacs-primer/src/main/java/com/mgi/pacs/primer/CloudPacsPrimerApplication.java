package com.mgi.pacs.primer;


import org.springframework.boot.SpringApplication;
import org.springframework.boot.autoconfigure.SpringBootApplication;
import springfox.documentation.swagger2.annotations.EnableSwagger2;

@EnableSwagger2
@SpringBootApplication
public class CloudPacsPrimerApplication {

    public static void main(String[] args) {

        SpringApplication app = new SpringApplication(CloudPacsPrimerApplication.class);
        app.addInitializers(new com.mgi.pacs.adapter.listenler.PcrApplicationRunner());
        app.run(args);
    }

}
