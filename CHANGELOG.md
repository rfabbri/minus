# Changelog

MINUS: MInimial problem NUmerical continuation Solver 
Copyright (C) 2018-2025 Ricardo Fabbri

All notable changes to MiNuS will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## 20260315

### Added
- MILESTONE: Tutorial for adding your own problem fully working and tested (see tutorial/README.md)
- commands in cmd/ modularized (easier to maintain and understand)
- better documetation for evaluators HxHt and HxH
- LLOG directives to print information that is shown only when compiling in
  Debug mode, while documeting functionality
- AGENTS.md and GEMINI.md files for coding with AI

## 20250717 Commit caef5e51e48ee95f10f4d90e128b5e4770cfb345

### Added

A new flag to enable/disable prefiltering of degenerate input (default on).

For instance, this close to degenrate problem from the synthetic data may still
be solved since homotopy continuation is very robust:

synthdata 35 99 75 1815 2558 3880 0 1 | minus-chicago -i --prefilter_degeneracy=no

But you can detect this early-on and just discard this problem (say, in
RANSAC):

synthdata 35 99 75 1815 2558 3880 0 1 | minus-chicago -i --prefilter_degeneracy=yes

## 20260430 

### Added
Compile flag -flto for link time optimization. The benefit seems a couple ms
shavings on prelimuinary tests, but it is worth adding as it logically makes
sense. It may slowdown Release compile time but we don't care about that,
we want raw runtime speed.
