include .env

print:
	@echo "DATA_PATH: $(DATA_PATH)"

_create-volumes-dir: # 创建数据的目录
	mkdir -p $(DATA_PATH)/.data/mysql \
			 $(DATA_PATH)/.data/redis \
			 $(DATA_PATH)/.data/xpcr/files

LOCAL_SERVICES := mysql redis
SPRING_PCR:= cloud-pacs-primer

BUILDER := docker buildx bake

LocalData := docker-compose --env-file .env -p xpcr-local -f docker-compose.yaml -f docker-compose.localdata.yaml

local_data_up: _create-volumes-dir
	$(LocalData) up -d $(LOCAL_SERVICES)
local_data_down:
	$(LocalData) down -v --remove-orphans

local_client:
	source scripts/export-env.sh .env;\
	source scripts/export-env.sh .env.local;\
	yarn start:client
local_server:
	source scripts/export-env.sh .env;\
    source scripts/export-env.sh .env.local;\
    yarn start:server

local_client_ubuntu:
	bash start_local_client.sh

local_server_ubuntu:
	bash start_local_server.sh

clean: ## clean and delete git ignore and dirty files
	git clean -fxd
install:
	yarn install

build-x-pcr:
	$(BUILDER) x-pcr

build-x-search:
	$(BUILDER) x-search

build-x-search-dev:
	$(BUILDER) x-search-dev

build-spring_pcr:
	mvn -f ./mgi-pcr-adapter/pom.xml install -Dmaven.skip.test=true;\
	mvn -f ./mgi-pcr-adapter/cloud-pacs-primer/pom.xml install -Dmaven.skip.test=true;

build: build-x-pcr build-x-search build-x-search-dev build-spring_pcr
	@echo "Build Done"

#pcr-adapter 
spring_pcr_start:
	mvn -f ./mgi-pcr-adapter/pom.xml install -Dmaven.skip.test=true;\
	mvn -f ./mgi-pcr-adapter/cloud-pacs-primer/pom.xml install -Dmaven.skip.test=true;\
	$(LocalData) up -d $(SPRING_PCR)

spring_pcr_remove:
	$(LocalData) stop cloud-pacs-primer;\
	$(LocalData) rm -fv cloud-pacs-primer;\
	docker rmi cloud-pacs-primer:dev

# show docker ps
docker-compose-ps:
	$(LocalData) ps


