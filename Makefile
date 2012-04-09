# ----------------------------------
# builds bamtools and executable
# ----------------------------------
all:
	cd bamtools && mkdir -p build && cd build && cmake .. && $(MAKE)
	cd src && $(MAKE)

clean:
	cd src && $(MAKE) clean
	rm -f bin/*

# ----------------------------------
# Phony Targets
# cite: http://linuxdevcenter.com/pub/a/linux/2002/01/31/make_intro.html?page=2
#  - "special rule .PHONY tells Make which targets are not files. This avoids
#     conflict with files of the same name, and improves performance"
# ----------------------------------
.PHONY: all clean