group "default" {
  targets = ["x-pcr", "x-search","x-search-dev"]
}


variable "IMAGE_REGISTRY" {
  default = "docker.io"
}

variable "SEMVER_FULL" {
  default = "v0.0.0-alpha"
}

variable "IMAGE_TAG" {
  default = "latest"
}

target "x-pcr" {
  context = "."
  dockerfile = "packaging/Dockerfile.x-pcr"
  args = {
    SEMVER_FULL = SEMVER_FULL
  }
  tags = ["${IMAGE_REGISTRY}/mgi-x/x-pcr:latest", "${IMAGE_REGISTRY}/mgi-x/x-pcr:${IMAGE_TAG}"]
}

target "x-search" {
  context = "."
  dockerfile = "packaging/Dockerfile.x-search"
  args = {
    SEMVER_FULL = SEMVER_FULL
  }
  tags = ["${IMAGE_REGISTRY}/mgi-x/x-search:latest", "${IMAGE_REGISTRY}/mgi-x/x-search:${IMAGE_TAG}"]
}

target "x-search-dev" {
  context = "."
  dockerfile = "packaging/Dockerfile.x-search-dev"
  args = {
    SEMVER_FULL = SEMVER_FULL
  }
  tags = ["${IMAGE_REGISTRY}/mgi-x/x-search-dev:latest", "${IMAGE_REGISTRY}/mgi-x/x-search-dev:${IMAGE_TAG}"]
}