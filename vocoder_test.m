% Austin Gilbert
% B EE 235 Spring 2020
% Final Project: N-Channel Vocoder
% vocoder_test.m
close all; clear all; clc;

inputfile = 'sent001.wav';
outputfile = 'output.wav';
bandnum = 6;
graphon = true;
soundon = true;
noiseon = true;

y = vocoder(inputfile, outputfile, bandnum, soundon, graphon, noiseon);