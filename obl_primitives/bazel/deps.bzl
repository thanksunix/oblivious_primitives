load("@bazel_tools//tools/build_defs/repo:http.bzl", "http_archive")

#--------------------------------------------------

def deps():
  _com_google_googletest()

def _com_google_googletest():
  if not native.existing_rule("com_google_googletest"):
    http_archive(
        name = "com_google_googletest",
        sha256 = "9bf1fe5182a604b4135edc1a425ae356c9ad15e9b23f9f12a02e80184c3a249c",
        type = "tar.gz",
        strip_prefix = "googletest-release-1.8.1",
        urls = [
            "https://github.com/google/googletest/archive/release-1.8.1.tar.gz",
        ],
    )
