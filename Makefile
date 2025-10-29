all : flash

TARGET:=stt
CH32FUN:=../../ch32fun/ch32fun
TARGET_MCU:=CH32V003

include $(CH32FUN)/ch32fun.mk
# include ../../ch32v003fun/ch32v003fun.mk

flash : cv_flash
clean : cv_clean
