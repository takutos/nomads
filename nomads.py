#!/usr/bin/python
"""
This file contains code for running Nomads simulation model: an agent-based model simulating pastoralist spatial movements.

Copyright 2015 Takuto Sakamoto (takutos@mac.com)
"""

import sys
import os
import csv
import math
from datetime import datetime

import matplotlib
matplotlib.use('TkAgg')

import pylab as PL
import random as RD
import scipy as SP
import gdal
from gdalconst import *

import pymas

RD.seed()
randState = RD.getstate()

monthsInYear = 12		# length of a year
yearsInGeneration = 1		# length of a generation
generationsInRun = 20000	# length of a simulation run
repetition = 1		# number of repetition of simulation runs

envName = 'NE Nigeria'	# name of the simulated drylands
dataFolder = 'NENGA'	# name of input data folder
spatiotemporalInput = 'VI'		# name of spatiotemporal input data file (vegetation index)
stInputYear = [str(year) for year in range(2005, 2015)]		# list of years of spatiotemporal input data
stInputBand = 1		# band of spatiotemporal input data to read
stAdjustment = True		# an option for adjusting spatiotemporal input data
spatialInput1 = 'LC'		# name of spatial input data file 1 (land classes)
sInput1Year = '2010'	# year of spatial input data 1
sInput1Band = 2		# band of spatial input data 1 to read
spatialInput2 = 'NoTsetse'		# name of spatial input data file 2 (distribution of tsetse flies)
#spatialInput2 = 'Morsitans'		# name of spatial input data file 2 (distribution of tsetse flies)
sInput2Year = '1999'	# year of spatial input data 2
sInput2Band = 1		# band of spatial input data 2 to read
stAdjustment = True		# an option for adjusting spatiotemporal input data
stAdjBaseClasses = [8, 9, 10]	# land classese whose vegetation indices are not to be adjusted
stAdjFactor = 0.5	# vegetation adjustment factor
s2AdjFactor = 1.0	# tsetse adjustment factor
carryingCapacity = 3	# number of agents which each site can support over a period of one year
degradationFactor = 1.0		# degree of land degradation that overpopulation causes (0.0 - 1.0)

population = 10		# population of Nomad agents
moveRange = 100		# distance that can be moved without incurring cost
grazeRange = 4		# range of neighborhood to be grazed
scoutRange = 100		# range of neighborhood to be scouted for information
scoutFrequency = 0.2	# scouting frequency in one month per agent
exchangeFrequency = 0.0		# frequency of information exchange in one month per agent
numOfAlts = 100		# maximum number of alternative routes to be considered in the adaptation phase
exploreTendency = 0.1	# tendency to allow major changes in the movement plan
movementCost = 100		# cost incurred when moving beyond a specified distance in one month
disruptionEffect = 10	# effect of a disruption event
adaptParams = [-0.5, -0.0, 0.000001, 0.01]		# adaptParams: parameters controlling the agent's update of the movement plan ([minimum potential normalized to 0.0, reference potential normalized to 1.0, adaptation noise for minimum potential, daptation noise for reference potential])
knowledgeDecay = 0.0	# frequency of losing each piece of knowledge during knowledge inheritance between generations

expropriation = 0.0		# rate of the lands to be expropriated
exprCategory = [12]		# list of rangeland classes to be targeted
exprExemption = []		# list of months during which the land access remains open
divRows = 0		# number of rangeland division in vertical direction
divCols = 0		# number of rangeland division in horizontal direction
resourceShare = True	# an option for introducing resource sharing among agents (if False, a randomly chosen agent gets exclusive resource access)
instSlack = 0.0		# institutional slack in implementing group/private access

display = True		# an option for displaying the landscape and agents
drawingMode = 2		# landscape information to be shown (0: LandClass, 1: Tsetse, 2: Resources, 3: Availability, 4: GrazeCount)
showRoute = True	# an option for showing the movement paths of agents
showRegime = False	# an option for showing the rangeland regime condition
outInterval = 120		# time interval for displaying the output

folderName = '_Simulation'		# name of an output data folder
settingOut = True	# an option for writing out information on simulation setting
processOut = True	# an option for writing out recorded data on the simulation process in a csv format
snapshot = True		# an option for capturing a snapshot of the situation where the average yearly resources obtained by agents is highest
snapshotFinal = True		# an option for capturing a snapshot of the final state of each run
succesiveShots = True	# an option for taking successive snapshots of a simulation run at a specified interval of generations
shotInterval = 100	# intervals (generations) at which successive snapshots are taken
countMap = True		# an option for writing out grazing count maps in a GeoTIFF format
#countStart = 1001	# generation at which grazing count begins

#set control variables
conVars = [generationsInRun, population, drawingMode, outInterval]
conSetting = [[generationsInRun, 'scale', 'int', 'Generations', 1, 20000, 1],
			[population, 'scale', 'int', 'Population', 1, 100, 1],
			[drawingMode, 'scale', 'int', 'Display Mode (0: LC 1: TT 2: RES 3: AVL 4: GC)', 0, 4, 1],
			[outInterval, 'scale', 'int', 'Display Interval (months)', 1, 1200, 1]]

class Time(object):
	"""Represents a parent class for managing the simulation time."""

class Env(object):
	"""Represents a parent class for controlling the simulation environment."""

class Agt(object):
	"""Represents a parent class for implementing agents."""

class Inst(object):
	"""Represents an institution controlling interactions among agents."""


class NomadTime(Time):
	"""Controls the time structure and change in the model."""

	def __init__(self, lengthYear=12, lengthGeneration=10, lengthRun=1000):
		"""Constructor of the time agent.

		Args:
			lengthYear: number of months in one year
			lengthGeneration: number of years in one generation
			lengthRun: number of generations in one simulation run
		"""
		self.monthsInYear = lengthYear
		self.yearsInGeneration = lengthGeneration
		self.generationsInRun = lengthRun
		self.month = 0
		self.year = 0
		self.generation = 0

	def Update(self):
		"""Updates each component (month, year, generation) of the time structure.

		Returns:
			False when the last generation has been reached, and the simulation should end.
		"""
		self.month += 1
		if self.monthsInYear == 1 or self.month % self.monthsInYear == 1:
			self.month = 1
			self.year += 1
			if self.yearsInGeneration == 1 or self.year % self.yearsInGeneration == 1:
				self.year = 1
				self.generation += 1
				if self.generation > self.generationsInRun:
					return False
		return True

	def TimeInMonth(self):
		"""Measures the current time in the total number of months passed from the start of simulation.

		Returns:
			Time in month.
		"""
		return self.month + (self.year - 1) * self.monthsInYear + (self.generation - 1) * self.yearsInGeneration * self.monthsInYear


class Drylands(Env):
	"""Represents a dryland landscape."""

	def __init__(self, name='NE_Nigeria', capacity=3, degradeFactor=1.0):
		"""Constructor of the environment.

		Args:
			name: name of the environment
			capacity: number of agents which each site can support over a period of one year
			degradeFactor: degree of land degradation that overpopulation causes (0.0 - 1.0)
		"""
		self.name = name
		self.carryingCapacity = capacity
		self.degradationFactor = degradeFactor

	def ReadData(self, data_dir=os.path.join('Input', 'NENGA'), spatiotemporalInput='VI', stInputYear=['2014'], stInputBand=1, spatialInput1='LC', sInput1Year='2012', sInput1Band=1, spatialInput2='Marsitans', sInput2Year='1999', sInput2Band=1, stAdjustment=False, stAdjBaseClasses=[8, 9, 10], stAdjFactor=0.5, s2AdjFactor=1.0, time=None):
		"""Reads spatial data files, builds a spatial makeup of the landscape.

		Args:
			data_dir: string directory name
			spatiotemporalInput: name of the spatiotemporal data (vegetation index) to read
			stInputYear: list of years of spatiotemporal input data
			stInputBand: band of spatiotemporal input data to read
			spatialInput1: name of the first spatial (temporally invariant) data (land classes) to read
			sInput1Year: year of the first spatial input data
			sInput1Band: band of the first spatial input data to read
			spatialInput2: name of the second spatial (temporally invariant) data (tsetse flies) to read
			sInput2Year: year of the second spatial input data
			sInput2Band: band of the second spatial input data to read
			stAdjustment: an option for adjusting spatiotemporal input data
			stAdjBaseClasses: land classese whose vegetation indices are not to be adjusted
			stAdjFactor: vegetation adjustment factor
			s2AdjFactor: tsetse adjustment factor
			time: a NomadTime agent
		"""
		self.period = len(stInputYear) * 12
		veg = []
		for year in stInputYear:
			for i in xrange(12):
				month = str(i + 1)
				if len(month) == 1:
					month = '0' + month
				filename = spatiotemporalInput + '_' + year + month + '.tif'
				filepath = os.path.join(data_dir, filename)
				dataset = gdal.Open(filepath, GA_ReadOnly)
				if year == stInputYear[0] and i == 0:
					self.width = dataset.RasterXSize
					self.height = dataset.RasterYSize
					self.projection = dataset.GetProjection()
					self.geotransform = dataset.GetGeoTransform()
				band = dataset.GetRasterBand(stInputBand)
				veg.append(SP.flipud(band.ReadAsArray()) * 0.0001)
		self.vegetation = SP.array(veg)
		filename = spatialInput1 + '_' + sInput1Year + '.tif'
		filepath = os.path.join(data_dir, filename)
		dataset = gdal.Open(filepath, GA_ReadOnly)
		band = dataset.GetRasterBand(sInput1Band)
		self.landClass = SP.flipud(band.ReadAsArray())
		filename = spatialInput2 + '_' + sInput2Year + '.tif'
		filepath = os.path.join(data_dir, filename)
		dataset = gdal.Open(filepath, GA_ReadOnly)
		band = dataset.GetRasterBand(sInput2Band)
		self.tsetse = SP.flipud(band.ReadAsArray()) * s2AdjFactor
		self.resources = self.vegetation[0]
		self.availability = SP.ones([self.height, self.width])
		self.grazedSites = SP.zeros([self.height, self.width], dtype=SP.int16)
		self.grazingCounts = SP.zeros([time.monthsInYear, self.height, self.width], dtype=SP.int16)
		if stAdjustment:
			adjustFilter = SP.ones([self.height, self.width], dtype=SP.bool8)
			for i in stAdjBaseClasses:
				adjustFilter = SP.logical_and(adjustFilter, self.landClass != i)
			for i in xrange(self.period):
				SP.putmask(self.vegetation[i], self.vegetation[i] < 0, self.vegetation[i] * 0)
				SP.putmask(self.vegetation[i], adjustFilter, self.vegetation[i] * stAdjFactor)

	def Update(self, time):
		"""Updates the landscape conditions.

		Args:
			time: a NomadTime agent
		"""
		self.resources = self.vegetation[time.TimeInMonth() % self.period - 1]
		if time.month == 1:
			self.grazedSites = SP.zeros([self.height, self.width], dtype=SP.int16)
			if time.year == 1:
				self.availability = SP.ones([self.height, self.width])
		SP.putmask(self.availability, self.grazedSites > self.carryingCapacity, SP.ones([self.height, self.width]) * (1.0 - self.degradationFactor))

	def Draw(self, mode='LandClass', time=None):
		"""Displays the landscape condition based on a selected variable.

		Args:
			mode: the landscape variable to be shown (0: LandCover, 1: Tsetse, 2: Resources, 3: Availability, 4: GrazeCount)
		Returns
			the displayed profile
			time: a NomadTime agent
		"""
		if mode == 'LandClass' or mode == 0:
			cond = 'LandClass'
			PL.pcolormesh(self.landClass, vmin = -2, vmax = 16, cmap = PL.cm.gist_earth)
		elif mode == 'Tsetse' or mode == 1:
			cond = 'Tsetse'
			PL.pcolormesh(self.tsetse, vmin = 0.0, vmax = 0.1, cmap = PL.cm.autumn_r)
		elif mode == 'Resources' or mode == 2:
			cond = 'Resources'
			PL.pcolormesh(self.resources, vmin = 0.0, vmax = 1.0, cmap = PL.cm.Greens)
		elif mode == 'Availability' or mode == 3:
			cond = 'Availability'
			PL.pcolormesh(self.availability, vmin = 0.0, vmax = 1.0, cmap = PL.cm.BrBG)
		elif mode == 'GrazeCount' or mode == 4:
			cond = 'GrazeCount'
			PL.pcolormesh(self.grazingCounts[time.month - 1], vmin = 0.0, vmax = time.yearsInGeneration * time.generation, cmap = PL.cm.autumn_r)
		return cond


class Nomad(Agt):
	"""Represents a nomadic pastoral unit moving with livestock over the landscape."""

	def __init__(self, id=0, lengthRoute=12, adaptInterval=10, mRange=5, gRange=0, sRange=5, sFreq=0.1, eFreq=0.1, kDecay=0.01, numAlts=100, exploration=0.1, env=None):
		"""Constructor of the agent.

		Args:
			id: a unique ID of the agent
			lengthRoute: length of the agent's movement plan
			adaptInterval: number of movement cycles before the route update occurs
			mRange: distance that can be moved without incurring costs
			gRange: range of neighborhood to be grazed
			sRange: range of neighborhood to be scouted for information
			sFreq: scouting frequency in one month per agent
			kDecay: frequency of losing each piece of knowledge during knowledge inheritance between generations
			eFreq: frequency of information exchange in one month per agent
			numAlts: maximum number of alternative routes to be considered in the adaptation phase
			exploration: tendency to allow major changes in the movement plan
			env: a landscape agent (Env object) to be linked with
		"""
		self.id = id
		self.lengthOfRoute = lengthRoute
		self.adaptationInterval = adaptInterval
		self.moveRange = mRange
		self.grazeRange = gRange
		self.scoutRange = sRange
		self.scoutFrequency = sFreq
		self.knowledgeDecay = kDecay
		self.exchangeFrequency = eFreq
		self.numOfAlts = numAlts
		self.exploreTendency = exploration
		self.env = env
		initX = RD.randint(0, env.width - 1)
		initY = RD.randint(0, env.height - 1)
		xs = [initX for i in xrange(self.lengthOfRoute)]
		ys = [initY for i in xrange(self.lengthOfRoute)]
		#xs = [RD.randint(0, env.width - 1) for i in xrange(self.lengthOfRoute)]
		#ys = [RD.randint(0, env.height - 1) for i in xrange(self.lengthOfRoute)]
		self.route = SP.array([xs, ys])
		self.x = self.route[0][0]
		self.y = self.route[1][0]
		self.grazingSites = []
		self.clock = 0
		self.cycle = 0
		self.generation = 0
		self.institution = []
		self.group = 0
		self.knowledge = [[] for i in xrange(self.lengthOfRoute)]
		self.resources = SP.zeros(self.lengthOfRoute)
		self.costs = SP.zeros(self.lengthOfRoute)
		self.disruptions = SP.zeros(self.lengthOfRoute)
		self.potential = 0.0
		if id == 0:
			self.color = (1.0, 0.0, 0.0)
		else:
			self.color = (RD.random(), RD.random(), RD.random())

	def Move(self, env):
		"""Moves to a monthly camping site according to the movement plan

		Args:
			env: a landscape agent (Env object) to be linked with
		"""
		self.clock += 1
		if self.lengthOfRoute == 1 or self.clock % self.lengthOfRoute == 1:
			self.clock = 1
			self.cycle += 1
			if self.adaptationInterval == 1 or self.cycle % self.adaptationInterval == 1:
				self.cycle = 1
				self.generation += 1
		self.x = self.route[0][self.clock - 1]
		self.y = self.route[1][self.clock - 1]
		self.grazingSites = []
		for i in xrange(-self.grazeRange, self.grazeRange + 1):
			curX = self.x + i
			if curX >= 0 and curX < env.width:
				for j in xrange(-self.grazeRange, self.grazeRange + 1):
					curY = self.y + j
					if curY >= 0 and curY < env.height and ((curX - self.x) ** 2 + (curY - self.y) ** 2 <= self.grazeRange ** 2):
						self.grazingSites.append([curX, curY])
						env.grazedSites[curY][curX] += 1
						env.grazingCounts[self.clock - 1][curY][curX] += 1
						inst = self.institution[0]
						inst.nomadMap[curY][curX][self.id] = 1

	def Update(self, moveCost=10, disruptEffect=10, adaptParams=[-0.5, -0.0, 0.000001, 0.01], env=None, noms=None, rec=None):
		"""Iterates the agent's update rule.

		Args:
			moveCost: cost incurred when moving faster than at a specified speed
			disruptEffect: effect of a disruption event
			adaptParams: : parameters controlling the agent's update of the movement plan ([minimum potential normalized to 0.0, reference potential normalized to 1.0, adaptation noise for minimum potential, daptation noise for reference potential])
			env: a landscape agent (Env object) to be linked with
			noms: list of Nomad agents
			rec: a recording agent (Recorder object)
		"""
		self.Graze(env=env)
		self.Scout(env=env)
		self.Exchange(noms=noms)
		if self.clock == self.lengthOfRoute and self.cycle == self.adaptationInterval:
			nextRoute, data = self.Adapt(moveCost=moveCost, disruptEffect=disruptEffect, adaptParams=adaptParams)
			self.Report(rec=rec)
			self.Renew(nextRoute, data)

	def Graze(self, env):
		"""Herding the livestock around the camping site for grazing resources.

		Args:
			env: a landscape agent (Env object) to be linked with
		"""
		res = 0.0
		disruption = 0.0
		for site in self.grazingSites:
			curX, curY = site[0], site[1]
			inst = self.institution[0]
			access = inst.DetermineAccess(agtID=self.id, month=self.clock, x=curX, y=curY, env=env, slack=instSlack)
			res += access * env.availability[curY][curX] * env.resources[curY][curX]
			if RD.random() < env.tsetse[curY][curX]:
				disruption = 1.0
		self.resources[self.clock - 1] = (self.resources[self.clock - 1] * (self.cycle - 1) + res / len(self.grazingSites)) / self.cycle
		self.disruptions[self.clock - 1] = (self.disruptions[self.clock - 1] * (self.cycle - 1) + disruption) / self.cycle
		if self.cycle == 1:
			if self.clock == 1:
				dist = (self.route[0][self.clock - 1] - self.route[0][self.lengthOfRoute - 1]) ** 2 + (self.route[1][self.clock - 1] - self.route[1][self.lengthOfRoute - 1]) ** 2
			else:
				dist = (self.route[0][self.clock - 1] - self.route[0][self.clock - 2]) ** 2 + (self.route[1][self.clock - 1] - self.route[1][self.clock - 2]) ** 2
			if dist > self.moveRange ** 2:
				self.costs[self.clock - 1] = 1.0

	def Scout(self, env):
		"""Stochastically scouts a nearby site for information on grazing resources.

		Args:
			env: a landscape agent (Env object) to be linked with
		"""
		if RD.random() < self.scoutFrequency:
			while 1:
				dX = RD.choice(range(-self.scoutRange, self.scoutRange + 1))
				dY = RD.choice(range(-self.scoutRange, self.scoutRange + 1))
				if dX == 0 and dY == 0:
					continue
				if dX ** 2 + dY ** 2 > self.scoutRange ** 2:
					continue
				tgtX = self.x + dX
				tgtY = self.y + dY
				if (tgtX >= 0 and tgtX < env.width) and (tgtY >= 0 and tgtY < env.height):
					break
			res = 0.0
			disruption = 0.0
			counter = 0
			for i in xrange(-self.grazeRange, self.grazeRange + 1):
				curX = tgtX + i
				if curX >= 0 and curX < env.width:
					for j in xrange(-self.grazeRange, self.grazeRange + 1):
						curY = tgtY + j
						if curY >= 0 and curY < env.height and ((curX - tgtX) ** 2 + (curY - tgtY) ** 2 <= self.grazeRange ** 2):
							inst = self.institution[0]
							access = inst.DetermineAccess(agtID=self.id, month=self.clock, x=curX, y=curY, env=env, slack=instSlack)
							res += access * env.availability[curY][curX] * env.resources[curY][curX]
							if RD.random() < env.tsetse[curY][curX]:
								disruption = 1.0
							counter += 1
			search = False
			for info in self.knowledge[self.clock - 1]:
				if info[0] == tgtX and info[1] == tgtY:
					search = True
					info[2] += 1
					info[3] = (info[3] * (info[2] - 1) + res / counter) / info[2]
					info[4] = (info[4] * (info[2] - 1) + disruption) / info[2]
					break
			if not search:
				info = [tgtX, tgtY, 1, res / counter, disruption]
				self.knowledge[self.clock - 1].append(info)

	def Exchange(self, noms):
		"""Stochastically visits a nearby agent for exchanging information on grazing resources.

		Args:
			noms: list of Nomad agents
		"""
		if RD.random() < self.exchangeFrequency and len(noms) > 1:
			for nom in RD.sample(noms, len(noms)):
				if nom.id != self.id and ((nom.x - self.x) ** 2 + (nom.y - self.y) ** 2 <= self.scoutRange ** 2):
					month1 = RD.randint(1, self.lengthOfRoute)
					month2 = RD.randint(1, nom.lengthOfRoute)
					if len(self.knowledge[month1 - 1]) > 0 and len(nom.knowledge[month2 - 1]) > 0:
						info1 = RD.choice(self.knowledge[month1 - 1])
						info2 = RD.choice(nom.knowledge[month2 - 1])
						search = False
						for info in nom.knowledge[month1 - 1]:
							if info[0] == info1[0] and info[1] == info1[1]:
								search = True
								info[2] += info1[2]
								info[3] = (info[3] * (info[2] - info1[2]) + info1[3] * info1[2]) / info[2]
								info[4] = (info[4] * (info[2] - info1[2]) + info1[4] * info1[2]) / info[2]
								break
						if not search:
							nom.knowledge[month1 - 1].append(info1)
						search = False
						for info in self.knowledge[month2 - 1]:
							if info[0] == info2[0] and info[1] == info2[1]:
								search = True
								info[2] += info2[2]
								info[3] = (info[3] * (info[2] - info2[2]) + info2[3] * info2[2]) / info[2]
								info[4] = (info[4] * (info[2] - info2[2]) + info2[4] * info2[2]) / info[2]
								break
						if not search:
							self.knowledge[month2 - 1].append(info2)
					break

	def Adapt(self, moveCost, disruptEffect, adaptParams):
		"""Generates the next movement plan in an adaptive manner.

		Args:
			moveCost: cost incurred when moving faster than at a specified speed
			disruptEffect: effect of a disruption event
			adaptParams: : parameters controlling the agent's update of the movement plan ([minimum potential normalized to 0.0, reference potential normalized to 1.0, adaptation noise for minimum potential, daptation noise for reference potential])
		Returns:
			nextRoute: the chosen movement plan for the next generation
			data: integretaed knowledge obtained in the current generation
		"""

		# integrate scouting knowledge with information obtained through the regular movement
		data = [[] for i in xrange(self.lengthOfRoute)]
		num = 1
		for i in xrange(self.lengthOfRoute):
			data[i].append([self.route[0][i], self.route[1][i], self.adaptationInterval, self.resources[i], self.disruptions[i]])
			for info in self.knowledge[i]:
				if info[0] == self.route[0][i] and info[1] == self.route[1][i]:
					data[i][0][3] = (data[i][0][2] * data[i][0][3] + info[2] * info[3]) / (data[i][0][2] + info[2])
					data[i][0][4] = (data[i][0][2] * data[i][0][4] + info[2] * info[4]) / (data[i][0][2] + info[2])
					data[i][0][2] += info[2]
				else:
					data[i].append(info)
			if num < self.numOfAlts:
				num *= len(data[i])
				if num > self.numOfAlts:
					num = self.numOfAlts

		# stochastically generate and evaluate alternative movement plans
		routes = [[0 for i in xrange(self.lengthOfRoute)]]
		potentials = SP.zeros(num)
		#self.potential = (-1 * SP.mean(self.resources) + moveCost * SP.mean(self.costs) + disruptEffect * SP.mean(self.disruptions) - adaptParams[0]) / (adaptParams[1] - adaptParams[0])
		#potentials[0] = self.potential
		potential = moveCost * SP.mean(self.costs)
		for i in xrange(self.lengthOfRoute):
			info = data[i][0]
			potential -= info[3] / self.lengthOfRoute
			potential += disruptEffect * info[4] / self.lengthOfRoute
		potentials[0] = (potential - adaptParams[0]) / (adaptParams[1] - adaptParams[0])
		self.potential = potentials[0]
		for i in xrange(num - 1):
			while 1:
				alt = []
				for j in xrange(self.lengthOfRoute):
					if len(data[j]) > 1 and RD.random() < self.exploreTendency:
						alt.append(RD.choice(range(1, len(data[j]))))
					else:
						alt.append(0)
				if routes.count(alt) == 0:
					routes.append(alt)
					break
			potential = 0
			for j in xrange(self.lengthOfRoute):
				info = data[j][alt[j]]
				potential -= info[3] / self.lengthOfRoute
				if j == 0:
					dist = (info[0] - data[self.lengthOfRoute - 1][alt[self.lengthOfRoute - 1]][0]) ** 2 + (info[1] - data[self.lengthOfRoute - 1][alt[self.lengthOfRoute - 1]][1]) ** 2
				else:
					dist = (info[0] - data[j - 1][alt[j - 1]][0]) ** 2 + (info[1] - data[j - 1][alt[j - 1]][1]) ** 2
				if dist > self.moveRange ** 2:
					potential += moveCost / self.lengthOfRoute
				potential += disruptEffect * info[4] / self.lengthOfRoute
			potentials[i + 1] = (potential - adaptParams[0]) / (adaptParams[1] - adaptParams[0])
		temperature = (adaptParams[3] - adaptParams[2]) * self.potential + adaptParams[2]
		probs = SP.exp(-(1 / temperature) * potentials)

		# stochastically select the next movement plan among the alternatives
		lottery = sum(probs) * RD.random()
		for k, prob in enumerate(probs):
			lottery -= prob
			if lottery <= 0:
				break
		alt = routes[k]
		xs = [0 for i in xrange(self.lengthOfRoute)]
		ys = [0 for i in xrange(self.lengthOfRoute)]
		nextRoute = SP.array([xs, ys])
		for i in xrange(self.lengthOfRoute):
			nextRoute[0][i] = data[i][alt[i]][0]
			nextRoute[1][i] = data[i][alt[i]][1]

		return nextRoute, data

	def Renew(self, nextRoute=None, data=None):
		"""Renew the state of the nomad by initializng some of its variables.

		Args:
			nextRoute: a new movement plan
			data: knowledge to be inherited
		"""
		if nextRoute is None:
			initX = RD.randint(0, env.width - 1)
			initY = RD.randint(0, env.height - 1)
			xs = [initX for i in xrange(self.lengthOfRoute)]
			ys = [initY for i in xrange(self.lengthOfRoute)]
			#xs = [RD.randint(0, env.width - 1) for i in xrange(self.lengthOfRoute)]
			#ys = [RD.randint(0, env.height - 1) for i in xrange(self.lengthOfRoute)]
			self.route = SP.array([xs, ys])
		else:
			self.route = nextRoute
		self.x = self.route[0][0]
		self.y = self.route[1][0]
		if data == None or self.knowledgeDecay >= 1.0:
			self.knowledge = [[] for i in xrange(self.lengthOfRoute)]
		else:
			if self.knowledgeDecay > 0.0:
				for i in xrange(self.lengthOfRoute):
					for info in data[i]:
						if RD.random() < self.knowledgeDecay:
							data[i].remove(info)
			self.knowledge = data
		self.resources = SP.zeros(self.lengthOfRoute)
		self.costs = SP.zeros(self.lengthOfRoute)
		self.disruptions = SP.zeros(self.lengthOfRoute)
		self.potential = 0

	def Report(self, rec):
		"""Reports the agent's state (routes, resources, costs, disruptions, potential) to the Recorder agent.

		Args:
			rec: a recording agent (Recorder object)
		"""
		rec.routes[self.id] = self.route
		rec.resources[self.id] = self.resources
		rec.costs[self.id] = self.costs
		rec.disruptions[self.id] = self.disruptions
		rec.potentials[self.id] = self.potential
		rec.colors[self.id] = self.color

	def Show(self, route=True, size=5):
		"""Displays the agent on 2D landscape map.

		Args:
			route: an option for showing the site locations described by the movement plan
			size: size of the agent icon
		"""
		if route:
			xs = [x + 0.5 for x in self.route[0]]
			xs.append(self.route[0][0] + 0.5)
			ys = [y + 0.5 for y in self.route[1]]
			ys.append(self.route[1][0] + 0.5)
			PL.plot(xs, ys, color=self.color, ls=':', marker='+', alpha=0.8)
			for i in xrange(self.lengthOfRoute):
				curX, curY = self.route[0][i], self.route[1][i]
				PL.annotate(s=str(i + 1), xy=(curX + 0.5, curY + 0.5), color=self.color, ha='center', va='center', alpha=0.8)
		PL.scatter([self.x + 0.5], [self.y + 0.5], s=size ** 2, color=self.color)


class RangelandRegime(Inst):
	"""Represents a rangeland institution governing the land access of Nomad agents."""

	def __init__(self, id, range, time, env, noms, resShare=False):
		"""Constructor of the institution.

		Args:
			id: a unique ID of the institution
			range: list of IDs of the agents complying the institution
			time: a NomadTime agent
			env: a landscape agent (an Env object)
			noms: list of Nomad agents
			resShare: an option for introducing resource sharing among agents
		"""
		self.id = id
		self.range = range
		self.resourceShare = resShare
		self.monthsInYear = time.monthsInYear
		self.width = env.width
		self.height = env.height
		self.landAccess = SP.ones([self.monthsInYear, self.height, self.width], dtype=SP.int16) * 999
		self.expropriatedSites = []
		self.numOfGroups = 0
		self.groupMembers = []
		for nom in noms:
			if nom.id in self.range:
				nom.institution.append(self)
		self.nomadMap = SP.zeros([self.height, self.width, len(noms)], dtype=SP.int8)

	def Update(self):
		"""Updates institutional variables."""
		self.nomadMap *= 0

	def DetermineAccess(self, agtID, month, x, y, env, slack=0.0):
		"""Determines a fraction of available resources to which the given agent have access.

		Args:
			agtID: ID of a Nomad agent
			month: month of a year under consideration
			x: x coordinate of the location under consideration
			y: y coordinate of the location under consideration
			env: a landscape agent (an Env object)
			slack: slack in implementing group/private access
		"""
		if agtID not in self.range:
			return 0.0
		if self.landAccess[month - 1][y][x] == -1:
			return 0.0
		if self.landAccess[month - 1][y][x] == 999:
			info = self.nomadMap[y][x]
			pop = SP.sum(info)
			if (info[agtID] == 1 and pop == 1) or (info[agtID] == 0 and pop == 0):
				return 1.0
			if self.resourceShare:
				if info[agtID] == 0:
					pop += 1
				return 1.0 / float(pop)
			else:
				IDs = [i for i in xrange(len(info)) if info[i] == 1]
				if agtID not in IDs:
					IDs.append(agtID)
				if agtID == RD.choice(IDs):
					return 1.0
				else:
					return 0.0
		k = self.landAccess[month - 1][y][x]
		info = SP.zeros(len(self.nomadMap[y][x]), SP.int8)
		for i in xrange(len(self.nomadMap[y][x])):
			if i in self.groupMembers[k - 1] and self.nomadMap[y][x][i] == 1:
				info[i] = 1
		pop = SP.sum(info)
		if agtID not in self.groupMembers[k - 1]:
			if RD.random() < slack:
				if self.resourceShare:
					return 1.0 / float(pop + 1)
				else:
					IDs = [i for i in xrange(len(info)) if info[i] == 1]
					IDs.append(agtID)
					if agtID == RD.choice(IDs):
						return 1.0
					else:
						return 0.0
			else:
				return 0.0
		else:
			if self.resourceShare:
				if info[agtID] == 0:
					pop += 1
				return 1.0 / float(pop)
			else:
				IDs = [i for i in xrange(len(info)) if info[i] == 1]
				if agtID not in IDs:
					IDs.append(agtID)
				if agtID == RD.choice(IDs):
					return 1.0
				else:
					return 0.0

	def ExpropriateLand(self, env, rate=0.1, category=None, exemption=None):
		"""Expropriates randomly chosen rangelands, and renders these areas inaccessible.

		Args:
			env: a landscape agent (an Env object)
			rate: rate of the lands to be expropriated
			category: list of integers denoting targeted land classes
			exemption: list of months during which the land access remains open
		"""
		if rate <= 0.0:
			return
		for i in xrange(self.width):
			for j in xrange(self.height):
				if category == None or category == [] or env.landClass[j][i] in category:
					if RD.random() < rate:
						for m in xrange(self.monthsInYear):
							if (m + 1) not in exemption:
								self.landAccess[m][j][i] = -1
								self.expropriatedSites.append([m + 1, i, j])
		self.expropriatedSites = SP.array(self.expropriatedSites)

	def DivideLand(self, noms, rows=0, columns=0):
		"""Devides the rangeland, and assings each tract to a group of Nomad agents.

		Args:
			noms: list of Nomad agents
			rows: number of division in vertical direction
			columns: number of division in horizontal direction
		"""
		if rows == 0 or columns == 0:
			return
		rInterval = float(self.height) / rows
		cInterval = float(self.width) / columns
		for m in xrange(self.monthsInYear):
			for i in xrange(self.width):
				for j in xrange(self.height):
					rNo = j // rInterval
					cNo = i // cInterval
					self.landAccess[m][j][i] = 1 + cNo + rNo * columns
		self.numOfGroups = rows * columns
		self.groupMembers = [[] for i in xrange(self.numOfGroups)]
		for nom in noms:
			nom.group = nom.id % self.numOfGroups + 1
			self.groupMembers[nom.group - 1].append(nom.id)

	def Show(self, time):
		"""Displays the rangeland regime (expropriated sites) on 2D landscape map.

		Args:
			time: a NomadTime agent
		"""
		for site in self.expropriatedSites:
			if site[0] == time.month:
				PL.plot(site[1] + 0.5, site[2] + 0.5, color=(1.0, 0.5, 0.0), marker='x', ms=5, alpha=0.5)


class Recorder(Agt):
	"""Collects and arranges various data on the simulation process, converting them to appropriate file formats."""

	def __init__(self, data_dir='Output', pop=1, lengthRoute=12, setting=True, process=True, snap=True, snapFin=True, shots=False, sInterval=100, cMap=True):
		"""Constructor of the agent.

		Args:
			data_dir: string directory name
			pop: population of Nomad agents
			lengthRoute: length of an agent's movement plan
			setting: if True, write out information on simulation setting
			process: if True, write out recorded data on the simulation process in a csv format
			snap: if True, record and capture a snapshot of the situation where the average yearly resources obtained by agents is highest.
			snapFin: if True, record and capture a snapshot of the final state of a simulation run.
			shots: taking successive snapshots of a simulation run at a specified interval of generations
			sInterval: intervals (generations) at which successive snapshots are taken
			cMap: if True, write out grazing count maps in a GeoTIFF format
		"""
		self.data_dir = data_dir
		self.population = pop
		self.lengthOfRoute = lengthRoute
		self.settingOut = setting
		self.processOut = process
		self.snapshot = snap
		self.snapshotFinal = snapFin
		self.succesiveShots = shots
		self.shotInterval = sInterval
		self.countMap = cMap
		self.routes = SP.zeros([self.population, 2, self.lengthOfRoute], dtype=SP.int16)
		self.resources = SP.zeros([self.population, self.lengthOfRoute])
		self.costs = SP.zeros([self.population, self.lengthOfRoute])
		self.disruptions = SP.zeros([self.population, self.lengthOfRoute])
		self.potentials = SP.zeros(self.population)
		self.colors = SP.zeros([self.population, 3])
		self.grossPotential = 0.0
		self.resourceTimeMeans = SP.zeros(self.population)
		self.resourceTimeStds = SP.zeros(self.population)
		self.resourceTimeAgtMean = 0.0
		self.resourceTimeMeanMax = 0.0
		self.resourceTimeMeanMedian = 0.0
		self.resourceTimeMeanMin = 0.0
		self.resourceTimeMeanStd = 0.0
		self.resourceTimeStdAgtMean = 0.0
		self.maxNomadID = -1
		self.medianNomadID = -1
		self.minNomadID = -1
		self.costTimeMeans = SP.zeros(self.population)
		self.costTimeAgtMean = 0.0
		self.disruptionTimeMeans = SP.zeros(self.population)
		self.disruptionTimeAgtMean = 0.0
		self.displacementVectors = SP.zeros([self.population, 2, self.lengthOfRoute])
		self.displacementAgtMeans = SP.zeros([2, self.lengthOfRoute])
		self.distances = SP.zeros([self.population, self.lengthOfRoute])
		self.distanceTimeMeans = SP.zeros(self.population)
		self.distanceAgtMeans = SP.zeros(self.lengthOfRoute)
		self.distanceTimeAgtMean = 0.0
		self.distanceTimeMeanStd = 0.0
		self.distanceMax = 0.0
		self.distanceMin = 0.0
		self.synchros = SP.zeros(self.lengthOfRoute)
		self.synchroTimeMean = 0.0
		self.centroids = SP.zeros([self.population, 2])
		self.centroidAgtMean = SP.array([0.0, 0.0])
		self.momentXs = SP.zeros(self.population)
		self.momentXAgtMean = 0.0
		self.momentXMax = 0.0
		self.momentXMin = 0.0
		self.momentXStd = 0.0
		self.momentYs = SP.zeros(self.population)
		self.momentYAgtMean = 0.0
		self.momentYMax = 0.0
		self.momentYMin = 0.0
		self.momentYStd = 0.0
		self.maxRanges = SP.zeros(self.population)
		self.maxRangeMean = 0.0
		self.maxRangeMax = 0.0
		self.maxRangeMin = 0.0
		self.maxRangeStd = 0.0
		self.maxNomadResources = 0.0
		self.maxNomadDisruptions = 0.0
		self.maxNomadDistance = 0.0
		self.maxNomadCentroid = SP.array([0.0, 0.0])
		self.maxNomadMomentX = 0.0
		self.maxNomadMomentY = 0.0
		self.maxNomadMaxRange = 0.0
		self.medianNomadResources = 0.0
		self.medianNomadDisruptions = 0.0
		self.medianNomadDistance = 0.0
		self.medianNomadCentroid = SP.array([0.0, 0.0])
		self.medianNomadMomentX = 0.0
		self.medianNomadMomentY = 0.0
		self.medianNomadMaxRange = 0.0
		self.minNomadResources = 0.0
		self.minNomadDisruptions = 0.0
		self.minNomadDistance = 0.0
		self.minNomadCentroid = SP.array([0.0, 0.0])
		self.minNomadMomentX = 0.0
		self.minNomadMomentY = 0.0
		self.minNomadMaxRange = 0.0
		self.landDegradation = 0.0
		self.recordedGeneration = -1
		self.recordedMeanResources = 0.0
		self.recordedRoutes = SP.zeros([self.population, 2, self.lengthOfRoute], dtype=SP.int16)
		self.recordedResources = SP.zeros(self.population)
		self.recordedCosts = SP.zeros(self.population)
		self.recordedDisruptions = SP.zeros(self.population)
		self.recordedVectors = SP.zeros([self.population, 2, self.lengthOfRoute])
		self.recordedMeanDistance = SP.zeros(self.population)
		self.recordedCentroids = SP.zeros([self.population, 2])
		self.recordedMomentXs = SP.zeros(self.population)
		self.recordedMomentYs = SP.zeros(self.population)
		self.recordedMaxRanges = SP.zeros(self.population)
		self.recordedLandDegradation = 0.0

	def Update(self, run=1, time=None, env=None, inst=None, mode='LandClass'):
		"""Arranges the simulation data and writes them out at the end of each generation.

		Args:
			run: number of the present simulation run
			time: a NomadTime agent
			env: a landscape agent (an Env object)
			inst: a land institution agent (a RangelandRegime object)
			mode: the landscape variable to be shown (0: LandClass, 1: Tsetse, 2: Resources, 3: Availability)
		"""
		if time.month < time.monthsInYear or time.year < time.yearsInGeneration:
			return

		# recording the setting specification
		if settingOut and run == 1 and time.generation == 1:
			filename = os.path.join(self.data_dir, 'setting.csv')
			fp = open(filename, 'a')
			dt = datetime.now()
			fp.write(dt.strftime("%Y%m%d%H%M") + '\n')
			line = 'EnvName,MonthsInYear,YearsInGeneration,GenerationsInRun,Repetition,SpatiotemporalInput,STInputYear,STInputBand,SpatialInput1,SInput1Year,SInput1Band,SpatialInput2,SInput2Year,SInput2Band,SpatiotemporalAdjustment,STAdjBaseClasses,STAdjFactor,S2AdjFactor,CarryingCapacity,DegradationFactor,Population,MoveRange,GrazeRange,ScoutRange,ScoutFrequency,ExchangeFrequency,NumOfAlts,ExploreTendency,MovementCost,DisruptionEffect,AdaptParams,KnowledgeDecay,Expropriation,ExprCategory,ExprExemption,DivRows,DivCols,ResourceShare,InstSlack\n'
			fp.write(line)
			line = envName + ',' + str(time.monthsInYear) + ',' + str(time.yearsInGeneration) + ',' + str(time.generationsInRun) + ',' + str(repetition) + ',' + str(spatiotemporalInput) + ',' + str(stInputYear).replace(' ', '').replace('\'', '').replace(',', ' ') + ',' + str(stInputBand) + ',' + str(spatialInput1) + ',' + str(sInput1Year) + ',' + str(sInput1Band) + ',' + str(spatialInput2) + ',' + str(sInput2Year) + ',' + str(sInput2Band) + ',' + str(stAdjustment) + ',' + str(stAdjBaseClasses).replace(' ', '').replace(',', ' ') + ',' + str(stAdjFactor) + ',' + str(s2AdjFactor) + ',' + str(carryingCapacity) + ',' + str(degradationFactor) + ',' + str(population) + ',' + str(moveRange) + ',' + str(grazeRange) + ',' + str(scoutRange) + ',' + str(scoutFrequency) + ',' + str(exchangeFrequency) + ',' + str(numOfAlts) + ',' + str(exploreTendency) + ',' + str(movementCost) + ',' + str(disruptionEffect) + ',' + '[' + str(adaptParams[0]) + ' ' + str(adaptParams[1]) + ' ' + str(adaptParams[2]) + ' ' + str(adaptParams[3]) + ']' + ',' + str(knowledgeDecay) + ',' + str(expropriation) + ',' + str(exprCategory).replace(' ', '').replace(',', ' ') + ',' + str(exprExemption).replace(' ', '').replace(',', ' ') + ',' + str(divRows) + ',' + str(divCols) + ',' + str(resourceShare) + ',' + str(instSlack) + '\n'
			fp.write(line)
			fp.close()
			filename = os.path.join(self.data_dir, 'seed.txt')
			fp = open(filename, 'a')
			fp.write(dt.strftime("%Y%m%d%H%M") + '\n')
			fp.write(str(randState) + '\n')
			fp.close()

		# calculating and writing out the output variables
		if processOut:
			self.grossPotential = sum(self.potentials)
			self.resourceTimeMeans = SP.mean(self.resources, 1)
			self.resourceTimeStds = SP.std(self.resources, 1)
			self.resourceTimeAgtMean = SP.mean(self.resourceTimeMeans)
			self.resourceTimeMeanMax = SP.amax(self.resourceTimeMeans)
			if self.population % 2 == 1:
				self.resourceTimeMeanMedian = SP.median(self.resourceTimeMeans)
			else:
				self.resourceTimeMeanMedian = SP.sort(self.resourceTimeMeans)[self.population / 2 - 1]
			self.resourceTimeMeanMin = SP.amin(self.resourceTimeMeans)
			self.resourceTimeMeanStd = SP.std(self.resourceTimeMeans)
			self.resourceTimeStdAgtMean = SP.mean(self.resourceTimeStds)
			self.maxNomadID = SP.argmax(self.resourceTimeMeans)
			for k, x in enumerate(self.resourceTimeMeans):
				if x == self.resourceTimeMeanMedian:
					break
			self.medianNomadID = k
			self.minNomadID = SP.argmin(self.resourceTimeMeans)
			self.costTimeMeans = SP.mean(self.costs, 1)
			self.costTimeAgtMean = SP.mean(self.costTimeMeans)
			self.disruptionTimeMeans = SP.mean(self.disruptions, 1)
			self.disruptionTimeAgtMean = SP.mean(self.disruptionTimeMeans)
			for i in xrange(self.population):
				xs = self.routes[i][0]
				ys = self.routes[i][1]
				self.maxRanges[i] = 0
				for j in xrange(self.lengthOfRoute):
					if j == self.lengthOfRoute - 1:
						self.displacementVectors[i][0][j] = xs[0] - xs[j]
						self.displacementVectors[i][1][j] = ys[0] - ys[j]
					else:
						self.displacementVectors[i][0][j] = xs[j + 1] - xs[j]
						self.displacementVectors[i][1][j] = ys[j + 1] - ys[j]
						for k in xrange(j + 1, self.lengthOfRoute):
							if (xs[k] - xs[j]) ** 2 + (ys[k] - ys[j]) ** 2 > self.maxRanges[i] ** 2:
								self.maxRanges[i] = math.sqrt((xs[k] - xs[j]) ** 2 + (ys[k] - ys[j]) ** 2)
				self.distances[i] = SP.sqrt(self.displacementVectors[i][0] ** 2 + self.displacementVectors[i][1] ** 2)
				self.centroids[i][0] = SP.mean(xs)
				self.centroids[i][1] = SP.mean(ys)
				self.momentXs[i] = SP.std(xs)
				self.momentYs[i] = SP.std(ys)
			self.displacementAgtMeans = SP.mean(self.displacementVectors, 0)
			self.distanceTimeMeans = SP.mean(self.distances, 1)
			self.distanceAgtMeans = SP.mean(self.distances, 0)
			self.distanceTimeAgtMean = SP.mean(self.distanceTimeMeans)
			self.distanceTimeMeanStd = SP.std(self.distanceTimeMeans)
			self.distanceMax = SP.amax(self.distances)
			self.distanceMin = SP.amin(self.distances)
			self.synchros = SP.sqrt((self.population * self.displacementAgtMeans[0]) ** 2 + (self.population * self.displacementAgtMeans[1]) ** 2) / (self.population * self.distanceAgtMeans)
			synchrosTrimmed = SP.array([x for x in self.synchros if x >= 0.0 and x <= 1.0])
			self.synchroTimeMean = SP.mean(synchrosTrimmed)
			self.centroidAgtMean = SP.mean(self.centroids, 0)
			self.momentXAgtMean = SP.mean(self.momentXs)
			self.momentXMax = SP.amax(self.momentXs)
			self.momentXMin = SP.amin(self.momentXs)
			self.momentXStd = SP.std(self.momentXs)
			self.momentYAgtMean = SP.mean(self.momentYs)
			self.momentYMax = SP.amax(self.momentYs)
			self.momentYMin = SP.amin(self.momentYs)
			self.momentYStd = SP.std(self.momentYs)
			self.maxRangeMean = SP.mean(self.maxRanges)
			self.maxRangeMax = SP.amax(self.maxRanges)
			self.maxRangeMin = SP.amin(self.maxRanges)
			self.maxRangeStd = SP.std(self.maxRanges)
			self.maxNomadResources = self.resourceTimeMeans[self.maxNomadID]
			self.maxNomadDisruptions = self.disruptionTimeMeans[self.maxNomadID]
			self.maxNomadDistance = self.distanceTimeMeans[self.maxNomadID]
			self.maxNomadCentroid = self.centroids[self.maxNomadID]
			self.maxNomadMomentX = self.momentXs[self.maxNomadID]
			self.maxNomadMomentY = self.momentYs[self.maxNomadID]
			self.maxNomadMaxRange = self.maxRanges[self.maxNomadID]
			self.medianNomadResources = self.resourceTimeMeans[self.medianNomadID]
			self.medianNomadDisruptions = self.disruptionTimeMeans[self.medianNomadID]
			self.medianNomadDistance = self.distanceTimeMeans[self.medianNomadID]
			self.medianNomadCentroid = self.centroids[self.medianNomadID]
			self.medianNomadMomentX = self.momentXs[self.medianNomadID]
			self.medianNomadMomentY = self.momentYs[self.medianNomadID]
			self.medianNomadMaxRange = self.maxRanges[self.medianNomadID]
			self.minNomadResources = self.resourceTimeMeans[self.minNomadID]
			self.minNomadDisruptions = self.disruptionTimeMeans[self.minNomadID]
			self.minNomadDistance = self.distanceTimeMeans[self.minNomadID]
			self.minNomadCentroid = self.centroids[self.minNomadID]
			self.minNomadMomentX = self.momentXs[self.minNomadID]
			self.minNomadMomentY = self.momentYs[self.minNomadID]
			self.minNomadMaxRange = self.maxRanges[self.minNomadID]
			self.landDegradation = 1.0 - SP.mean(env.availability)
			filename = os.path.join(self.data_dir, 'dynamics.csv')
			fp = open(filename, 'a')
			if time.generation == 1:
				dt = datetime.now()
				line = 'Run' + str(run) + '(' + dt.strftime("%Y%m%d%H%M") + ')\n'
				fp.write(line)
				line = 'Generation,Year,Month,GrossPotential,ResourceTimeAgtMean,ResourceTimeMeanMax,ResourceTimeMeanMedian,ResourceTimeMeanMin,ResourceTimeMeanStd,ResourceTimeStdAgtMean,CostTimeAgtMean,DisruptionTimeAgtMean,DisplacementAgtMeansX,DisplacementAgtMeansY,DistanceAgtMeans,DistanceTimeAgtMean,DistanceTimeMeanStd,DistanceMax,DistanceMin,Synchros,SynchroTimeMean,CentroidAgtMean,MomentXAgtMean,MomentXMax,MomentXMin,MomentXStd,MomentYAgtMean,MomentYMax,MomentYMin,MomentYStd,MaxRangeMean,MaxRangeMax,MaxRangeMin,MaxRangeStd,MaxNomadID,MaxNomadResources,MaxNomadDisruptions,MaxNomadDistance,MaxNomadCentroid,MaxNomadMomentX,MaxNomadMomentY,MaxNomadMaxRange,MedianNomadID,MedianNomadResources,MedianNomadDisruptions,MedianNomadDistance,MedianNomadCentroid,MedianNomadMomentX,MedianNomadMomentY,MedianNomadMaxRange,MinNomadID,MinNomadResources,MinNomadDisruptions,MinNomadDistance,MinNomadCentroid,MinNomadMomentX,MinNomadMomentY,MinNomadMaxRange,LandDegradation\n'
				fp.write(line)
			line = str(time.generation) + ',' + str(time.year) + ',' + str(time.month) + ',' + str(self.grossPotential) + ',' + str(self.resourceTimeAgtMean) + ',' + str(self.resourceTimeMeanMax) + ',' + str(self.resourceTimeMeanMedian) + ',' + str(self.resourceTimeMeanMin) + ',' + str(self.resourceTimeMeanStd) + ',' + str(self.resourceTimeStdAgtMean) + ',' + str(self.costTimeAgtMean) + ',' + str(self.disruptionTimeAgtMean) + ',' + str(self.displacementAgtMeans[0]).replace(',', '').replace('[', '').replace(']', '').replace('\n', '').strip() + ',' + str(self.displacementAgtMeans[1]).replace(',', '').replace('[', '').replace(']', '').replace('\n', '').strip() + ',' + str(self.distanceAgtMeans).replace(',', '').replace('[', '').replace(']', '').replace('\n', '').strip() + ',' + str(self.distanceTimeAgtMean) + ',' + str(self.distanceTimeMeanStd) + ',' + str(self.distanceMax) + ',' + str(self.distanceMin) + ',' + str(self.synchros).replace(',', '').replace('[', '').replace(']', '').replace('\n', '').strip() + ',' + str(self.synchroTimeMean) + ',' + str(self.centroidAgtMean).replace(',', '').replace('[', '').replace(']', '').replace('\n', '').strip() + ',' + str(self.momentXAgtMean) + ',' + str(self.momentXMax) + ',' + str(self.momentXMin) + ',' + str(self.momentXStd) + ',' + str(self.momentYAgtMean) + ',' + str(self.momentYMax) + ',' + str(self.momentYMin) + ',' + str(self.momentYStd) + ',' + str(self.maxRangeMean) + ',' + str(self.maxRangeMax) + ',' + str(self.maxRangeMin) + ',' + str(self.maxRangeStd) + ',' + str(self.maxNomadID) + ',' + str(self.maxNomadResources) + ',' + str(self.maxNomadDisruptions) + ',' + str(self.maxNomadDistance) + ',' + str(self.maxNomadCentroid).replace(',', '').replace('[', '').replace(']', '').replace('\n', '').strip() + ',' + str(self.maxNomadMomentX) + ',' + str(self.maxNomadMomentY) + ',' + str(self.maxNomadMaxRange) + ',' + str(self.medianNomadID) + ',' + str(self.medianNomadResources) + ',' + str(self.medianNomadDisruptions) + ',' + str(self.medianNomadDistance) + ',' + str(self.medianNomadCentroid).replace(',', '').replace('[', '').replace(']', '').replace('\n', '').strip() + ',' + str(self.medianNomadMomentX) + ',' + str(self.medianNomadMomentY) + ',' + str(self.medianNomadMaxRange) + ',' + str(self.minNomadID) + ',' + str(self.minNomadResources) + ',' + str(self.minNomadDisruptions) + ',' + str(self.minNomadDistance) + ',' + str(self.minNomadCentroid).replace(',', '').replace('[', '').replace(']', '').replace('\n', '').strip() + ',' + str(self.minNomadMomentX) + ',' + str(self.minNomadMomentY) + ',' + str(self.minNomadMaxRange) + ',' + str(self.landDegradation) + '\n'
			fp.write(line)
			fp.close()
			# possibly adding a procedure for writing out a record on each agent at each generation

		# recording the state of agents when the average yearly resources obtained hit the record
		if self.snapshot:
			if SP.mean(self.resources) >= self.recordedMeanResources:
				self.recordedGeneration = time.generation
				self.recordedMeanResources = SP.mean(self.resources)
				self.recordedRoutes = SP.array(self.routes)
				self.recordedResources = SP.mean(self.resources, 1)
				self.recordedCosts = SP.mean(self.costs, 1)
				self.recordedDisruptions = SP.mean(self.disruptions, 1)
				self.recordedLandDegradation = 1.0 - SP.mean(env.availability)
				if self.processOut:
					self.recordedVectors = SP.array(self.displacementVectors)
					self.recordedMeanDistance = SP.array(self.distanceTimeMeans)
					self.recordedCentroids = SP.array(self.centroids)
					self.recordedMomentXs = SP.array(self.momentXs)
					self.recordedMomentYs = SP.array(self.momentYs)
					self.recordedMaxRanges = SP.array(self.maxRanges)
			if time.generation == time.generationsInRun:
				dt = datetime.now()
				filename = os.path.join(self.data_dir, 'snapshotRecs.csv')
				fp = open(filename, 'a')
				if run == 1:
					line = 'Time,Run,Generation,MeanResources,'
					for i in xrange(self.population):
						nom = 'Nomad' + str(i)
						line = line + nom + 'Route,' + nom + 'Resources,' + nom + 'Cost,' + nom + 'Disruptions,' + nom + 'Vecotrs,' + nom + 'Distance,' + nom + 'Centroid,' + nom + 'MomentX,' + nom + 'MomentY,' + nom + 'MaxRange,'
					line = line + 'LandDegradation'
					fp.write(line + '\n')
				line = dt.strftime("%Y%m%d%H%M") + ',' + str(run) + ',' + str(self.recordedGeneration) + ',' + str(self.recordedMeanResources) + ','
				for i in xrange(self.population):
					line = line + str(zip(self.recordedRoutes[i][0], self.recordedRoutes[i][1])).replace(',', '').replace('[', '').replace(']', '').replace('\n', '').strip() + ',' + str(self.recordedResources[i]) + ',' + str(self.recordedCosts[i]) + ',' + str(self.recordedDisruptions[i]) + ',' + str(zip(self.recordedVectors[i][0], self.recordedVectors[i][1])).replace(',', '').replace('[', '').replace(']', '').replace('\n', '').strip() + ',' + str(self.recordedMeanDistance[i]) + ',' + str(self.recordedCentroids[i]).replace(',', '').replace('[', '').replace(']', '').replace('\n', '').strip() + ',' + str(self.recordedMomentXs[i]) + ',' + str(self.recordedMomentYs[i]) + ',' + str(self.recordedMaxRanges[i]) + ','
				line = line + str(self.recordedLandDegradation)
				fp.write(line + '\n')
				fp.close()
				PL.cla()
				cond = env.Draw(mode=drawingMode, time=time)
				PL.hold(True)
				for site in inst.expropriatedSites:
					if site[0] == time.month:
						PL.plot(site[1] + 0.5, site[2] + 0.5, color=(1.0, 0.5, 0.0), marker='x', ms=5, alpha=0.5)
				for i in xrange(self.population):
					xs = [x + 0.5 for x in self.recordedRoutes[i][0]]
					xs.append(self.recordedRoutes[i][0][0] + 0.5)
					ys = [y + 0.5 for y in self.recordedRoutes[i][1]]
					ys.append(self.recordedRoutes[i][1][0] + 0.5)
					PL.plot(xs, ys, color=self.colors[i], ls=':', marker='o', ms=10, mec=self.colors[i], mfc=self.colors[i], alpha=1.0)
					for j in xrange(self.lengthOfRoute):
						curX, curY = self.recordedRoutes[i][0][j], self.recordedRoutes[i][1][j]
						PL.annotate(s=str(j + 1), xy=(curX + 0.5, curY + 0.5), color='w', size='xx-small', ha='center', va='center', alpha=1.0)
				PL.hold(False)
				PL.axis('off')
				PL.xlim([0, env.width - 1])
				PL.ylim([0, env.height - 1])
				PL.axis('image')
				PL.tick_params(axis='both', length=0)
				filename = os.path.join(self.data_dir, '{name}({cond})_rec_run{run}_gen{gen}_{time}.png'.format(name=env.name, cond=cond, run=run, gen=self.recordedGeneration, time=dt.strftime("%Y%m%d%H%M")))
				PL.savefig(filename, format='png', frameon=False)

		# recording the state of agents in the final generation
		if self.snapshotFinal:
			if time.generation == time.generationsInRun:
				dt = datetime.now()
				filename = os.path.join(self.data_dir, 'snapshotFinRecs.csv')
				fp = open(filename, 'a')
				if run == 1:
					line = 'Time,Run,Generation,MeanResources,'
					for i in xrange(self.population):
						nom = 'Nomad' + str(i)
						line = line + nom + 'Route,' + nom + 'Resources,' + nom + 'Cost,' + nom + 'Disruptions,' + nom + 'Vecotrs,' + nom + 'Distance,' + nom + 'Centroid,' + nom + 'MomentX,' + nom + 'MomentY,' + nom + 'MaxRange,'
					line = line + 'LandDegradation'
					fp.write(line + '\n')
				line = dt.strftime("%Y%m%d%H%M") + ',' + str(run) + ',' + str(time.generation) + ',' + str(SP.mean(self.resources)) + ','
				for i in xrange(self.population):
					line = line + str(zip(self.routes[i][0], self.routes[i][1])).replace(',', '').replace('[', '').replace(']', '').replace('\n', '').strip() + ',' + str(SP.mean(self.resources[i])) + ',' + str(SP.mean(self.costs[i])) + ',' + str(SP.mean(self.disruptions[i])) + ',' + str(zip(self.displacementVectors[i][0], self.displacementVectors[i][1])).replace(',', '').replace('[', '').replace(']', '').replace('\n', '').strip() + ',' + str(self.distanceTimeMeans[i]) + ',' + str(self.centroids[i]).replace(',', '').replace('[', '').replace(']', '').replace('\n', '').strip() + ',' + str(self.momentXs[i]) + ',' + str(self.momentYs[i]) + ',' + str(self.maxRanges[i]) + ','
				line = line + str(1.0 - SP.mean(env.availability))
				fp.write(line + '\n')
				fp.close()
				PL.cla()
				cond = env.Draw(mode=drawingMode, time=time)
				PL.hold(True)
				for site in inst.expropriatedSites:
					if site[0] == time.month:
						PL.plot(site[1] + 0.5, site[2] + 0.5, color=(1.0, 0.5, 0.0), marker='x', ms=5, alpha=0.5)
				for i in xrange(self.population):
					xs = [x + 0.5 for x in self.routes[i][0]]
					xs.append(self.routes[i][0][0] + 0.5)
					ys = [y + 0.5 for y in self.routes[i][1]]
					ys.append(self.routes[i][1][0] + 0.5)
					PL.plot(xs, ys, color=self.colors[i], ls=':', marker='o', ms=10, mec=self.colors[i], mfc=self.colors[i], alpha=1.0)
					for j in xrange(self.lengthOfRoute):
						curX, curY = self.routes[i][0][j], self.routes[i][1][j]
						PL.annotate(s=str(j + 1), xy=(curX + 0.5, curY + 0.5), color='w', size='xx-small', ha='center', va='center', alpha=1.0)
				PL.hold(False)
				PL.axis('off')
				PL.xlim([0, env.width - 1])
				PL.ylim([0, env.height - 1])
				PL.axis('image')
				PL.tick_params(axis='both', length=0)
				filename = os.path.join(self.data_dir, '{name}({cond})_finRec_run{run}_gen{gen}_{time}.png'.format(name=env.name, cond=cond, run=run, gen=time.generation, time=dt.strftime("%Y%m%d%H%M")))
				PL.savefig(filename, format='png', frameon=False)

		# taking successive snapshots of a simulation run at a specified interval
		if self.succesiveShots:
			if time.generation == 1 or time.generation % self.shotInterval == 0:
				dt = datetime.now()
				PL.cla()
				cond = env.Draw(mode=drawingMode, time=time)
				PL.hold(True)
				for i in xrange(self.population):
					xs = [x + 0.5 for x in self.routes[i][0]]
					xs.append(self.routes[i][0][0] + 0.5)
					ys = [y + 0.5 for y in self.routes[i][1]]
					ys.append(self.routes[i][1][0] + 0.5)
					PL.plot(xs, ys, color=self.colors[i], ls=':', marker='o', ms=10, mec=self.colors[i], mfc=self.colors[i], alpha=1.0)
					for j in xrange(self.lengthOfRoute):
						curX, curY = self.routes[i][0][j], self.routes[i][1][j]
						PL.annotate(s=str(j + 1), xy=(curX + 0.5, curY + 0.5), color='w', size='xx-small', ha='center', va='center', alpha=1.0)
				PL.hold(False)
				PL.axis('off')
				PL.xlim([0, env.width - 1])
				PL.ylim([0, env.height - 1])
				PL.axis('image')
				PL.tick_params(axis='both', length=0)
				filename = os.path.join(self.data_dir, '{name}({cond})_run{run}_gen{gen}_{time}.png'.format(name=env.name, cond=cond, run=run, gen=time.generation, time=dt.strftime("%Y%m%d%H%M")))
				PL.savefig(filename, format='png', frameon=False)

		# writing out grazing count maps
		if self.countMap:
			if time.generation == time.generationsInRun:
				dt = datetime.now()
				filename = os.path.join(self.data_dir, '{name}(GrazeCount)_run{run}_gen{gen}_{time}.tif'.format(name=env.name, run=run, gen=time.generation, time=dt.strftime("%Y%m%d%H%M")))
				driver = gdal.GetDriverByName('GTiff')
				dataset = driver.Create(filename, env.width, env.height, time.monthsInYear + 1, gdal.GDT_Float32)
				dataset.SetGeoTransform(env.geotransform)
				dataset.SetProjection(env.projection)
				for i in xrange(time.monthsInYear):
					dataset.GetRasterBand(i + 1).WriteArray(SP.flipud(env.grazingCounts[i] / float(time.generationsInRun * time.yearsInGeneration)))
				dataset.GetRasterBand(time.monthsInYear + 1).WriteArray(SP.flipud(SP.sum(env.grazingCounts, 0) / float(time.generationsInRun * time.yearsInGeneration)))
				dataset = None

		print 'generation:', time.generation
		print 'gross potential  mean resources  mean cost  mean disruptions  land degradation'
		print sum(self.potentials), SP.mean(self.resources), SP.mean(self.costs), SP.mean(self.disruptions), 1.0 - SP.mean(env.availability)


def init():
	global run, time, env, noms, inst, rec
	global repetition, generationsInRun, population, drawingMode, outInterval
	generationsInRun, population, drawingMode, outInterval = tuple(conVars)

	if settingOut or processOut or snapshot or succesiveShots:
		path = os.path.join('Output', folderName)
		if rpt == repetition:
			if not os.path.exists(path):
				os.mkdir(path)
	else:
		path = 'Output'
	run = True
	time = NomadTime(monthsInYear, yearsInGeneration, generationsInRun)
	env = Drylands(name=envName, capacity=carryingCapacity, degradeFactor=degradationFactor)
	env.ReadData(data_dir=os.path.join('Input', dataFolder), spatiotemporalInput=spatiotemporalInput, stInputYear=stInputYear, stInputBand=stInputBand, spatialInput1=spatialInput1, sInput1Year=sInput1Year, sInput1Band=sInput1Band, spatialInput2=spatialInput2, sInput2Year=sInput2Year, sInput2Band=sInput2Band, stAdjustment=stAdjustment, stAdjBaseClasses=stAdjBaseClasses, stAdjFactor=stAdjFactor, s2AdjFactor=s2AdjFactor, time=time)
	noms = [Nomad(id=i, lengthRoute=monthsInYear, adaptInterval=yearsInGeneration, mRange=moveRange, gRange=grazeRange, sRange=scoutRange, sFreq=scoutFrequency, eFreq=exchangeFrequency, kDecay=knowledgeDecay, numAlts=numOfAlts, exploration=exploreTendency, env=env) for i in xrange(population)]
	inst = RangelandRegime(id=0, range=range(population), resShare=resourceShare, time=time, env=env, noms=noms)
	inst.DivideLand(rows=divRows, columns=divCols, noms=noms)
	inst.ExpropriateLand(env, rate=expropriation, category=exprCategory, exemption=exprExemption)
	rec = Recorder(data_dir=path, pop=population, lengthRoute=monthsInYear, setting=settingOut, process=processOut, snap=snapshot, snapFin=snapshotFinal, shots=succesiveShots, sInterval=shotInterval, cMap=countMap)


def draw():
	global rpt, run, time, env, noms, inst
	global generationsInRun, population, drawingMode, outInterval
	#dummy1, dummy2, dummy3, dummy4 = tuple(conVars)

	if display and run and (outInterval == 1 or time.month == 0 or time.TimeInMonth() % outInterval == 0):
		PL.cla()
		cond = env.Draw(mode=drawingMode, time=time)
		PL.hold(True)
		if showRegime:
			inst.Show(time=time)
		for nom in noms:
			nom.Show(route=showRoute)
		PL.hold(False)
		PL.axis('off')
		PL.xlim([0, env.width - 1])
		PL.ylim([0, env.height - 1])
		PL.axis('image')
		PL.tick_params(axis='both', length=0)
		PL.title('{name} ({cond}): run={run} gen={gen} year={year} month={month}'.format(name=env.name, cond=cond, run=repetition+1-rpt, gen=time.generation, year=time.year, month=time.month))


def step():
	global gui, rpt, run, time, env, noms, inst, rec
	global generationsInRun, population, drawingMode, outInterval
	#dummy1, dummy2, dummy3, dummy4 = tuple(conVars)

	if run:
		run = time.Update()
		if not run:
			rpt -= 1
			if rpt <= 0:
				gui.runEvent()
			else:
				gui.resetModel()
				gui.runEvent()
		env.Update(time)
		for nom in noms:
			nom.Move(env=env)
		for nom in noms:
			nom.Update(moveCost=movementCost, disruptEffect=disruptionEffect, adaptParams=adaptParams, env=env, noms=noms, rec=rec)
		inst.Update()
		rec.Update(run=repetition+1-rpt, time=time, env=env, inst=inst, mode=drawingMode)


def main(name, *args):
	global gui, rpt

	rpt = repetition
	gui = pymas.GUI(title='Nomads Model', conVars=conVars, conSetting=conSetting)
	gui.start(func=[init,draw,step])

if __name__ == '__main__':
	main(*sys.argv)
