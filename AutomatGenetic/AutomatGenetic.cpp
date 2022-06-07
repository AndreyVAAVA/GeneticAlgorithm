// AutomatGenetic.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <random>
#include <list>
#include <algorithm>

class Individ {
public:
	double start;
	double end;
	double x;
	double y;
	double score;
	long mutationSteps;

public:

	bool operator< (const Individ& ind) const {
		if (score - ind.score < 0) return true;
		else return false;
	}

	Individ(double start, double end, long mutationSteps) {
		this->start = start;
		this->end = end;
		this->mutationSteps = mutationSteps;
		this->score = 0;
		std::random_device rd; // obtain a random number from hardware
		std::mt19937 gen(rd()); // seed the generator
		std::uniform_real_distribution<> distr(start, end);
		this->x = distr(gen);
		this->y = distr(gen);
	}

	void calculateSelf() {
		score = function(x, y);
	}


	double function(double x, double y) {
		return x * x + y * y;
	}

	void mutate() {
		double delta = 0;
		std::random_device rd;
		std::default_random_engine eng(rd());
		std::uniform_real_distribution<float> distr(0, 1);
		for (int i = 1; i <= mutationSteps; i++)
		{
			if (distr(eng) < 1 / mutationSteps)
			{
				delta += 1 / (pow(2, i));
			}
		}
		if (rand() % 2 == 1)
		{
			delta = end * delta;
		}
		else
		{
			delta = start * delta;

		}
		x += delta;
		if (x < 0)
		{
			x = fmaxf(x, start);
		}
		else
		{
			x = fminf(x, end);
		}
		delta = 0;
		for (int j = 1; j < mutationSteps + 1; j++)
		{
			if (distr(eng) < 1 / mutationSteps)
			{
				delta += 1 / pow(2, j);
			}
		}
		if (rand() % 2 == 1)
		{
			delta = end * delta;
		}
		else
		{
			delta = start * delta;

		}
		y += delta;
		if (y < 0)
		{
			y = fmaxf(y, start);
		}
		else
		{
			y = fminf(y, end);
		}
	}
};

class Genetic
{
private:
	long numberOfIndividums;
	double crossoverRate;
	long mutationSteps;
	double chanceMutations;
	long numberLives;
	long start;
	long end;
	double bestScore;
	double* xy;
public:
	Genetic(long numberOfIndividums,
		double crossoverRate, long mutationSteps,
		double chanceMutations, long numberLives,
		long start, long end)
	{
		this->numberOfIndividums = numberOfIndividums;
		this->crossoverRate = crossoverRate;
		this->mutationSteps = mutationSteps;
		this->chanceMutations = chanceMutations;
		this->numberLives = numberLives;
		this->start = start;
		this->end = end;
		this->bestScore = std::numeric_limits<double>::infinity();
		this->xy = new double[2];
		xy[0] = std::numeric_limits<double>::infinity();
		xy[1] = std::numeric_limits<double>::infinity();
	}
	~Genetic()
	{

	}

	double getBestScore() {
		return bestScore;
	}

	double* getXY() {
		return xy;
	}

	std::vector<Individ> crossover(Individ parent1, Individ parent2)
	{
		Individ child1 = Individ(start, end, mutationSteps);
		Individ child2 = Individ(start, end, mutationSteps);
		std::random_device rand_dev;
		std::mt19937 generator(rand_dev());
		std::uniform_real_distribution<double> alpha = std::uniform_real_distribution<double>(0.01, 1);
		child1.x = parent1.x + alpha(generator) * (parent2.x - parent1.x);

		alpha = std::uniform_real_distribution<double>(0.01, 1);
		child1.y = parent1.y + alpha(generator) * (parent2.y - parent1.y);

		alpha = std::uniform_real_distribution<double>(0.01, 1);
		child2.x = parent1.x + alpha(generator) * (parent1.x - parent2.x);

		alpha = std::uniform_real_distribution<double>(0.01, 1);
		child2.y = parent1.y + alpha(generator) * (parent1.y - parent2.y);
		std::vector<Individ> vect;
		vect.push_back(child1);
		vect.push_back(child2);
		return vect;
	}

	void startGenetic()
	{
		std::vector<double> pack;
		pack.push_back(start);
		pack.push_back(end);
		pack.push_back(mutationSteps);
		std::vector<Individ> population;
		long i = 0;
		long j = 0;
		for (; i < numberOfIndividums; i++)
		{
			population.push_back(Individ(pack[0], pack[1], ceil(pack[2])));
		}
		srand(time(NULL));
		for (i = 0; i < numberLives; i++)
		{
			std::sort(population.begin(), population.end());
			std::vector<Individ> bestPopulation;
			std::vector<long> elem;
			for (; j < numberOfIndividums * crossoverRate; j++)
			{
				bestPopulation.push_back(population[j]);
				elem.push_back(j);
			}
			std::vector<Individ> child;
			for (Individ& individ1 : bestPopulation) {

				Individ& individ2 = bestPopulation[rand() % bestPopulation.size()];
				while (&individ1 == &individ2)
				{
					individ2 = bestPopulation[rand() % bestPopulation.size()];

				}
				std::vector<Individ> chd = crossover(individ1, individ2);
				child.push_back(chd[0]);
				child.push_back(chd[1]);
			}
			j = 0;
			for (j = 0; j < child.size(); j++)
			{
				population.push_back(child[j]);
			}
			for (Individ& individ : population) {
				individ.mutate();
				individ.calculateSelf();
			}
			std::sort(population.begin(), population.end());
			std::vector<Individ> nPopulation;
			for (j = 0; j < numberOfIndividums; j++)
			{
				nPopulation.push_back(population[j]);
			}
			population = nPopulation;
			if (population[0].score < bestScore)
			{
				bestScore = population[0].score;
				xy[0] = population[0].x;
				xy[1] = population[0].y;
			}
		}

	}
};

int main()
{
	setlocale(LC_ALL, "Ru");
	Genetic* a = new Genetic(500, 0.5, 15, 0.4, 200, -5, 5);
	a->startGenetic();
	std::cout << "ÎÏÒÈÌÈÇÈÐÎÂÀÍÍÎÅ ÇÍÀ×ÅÍÈÅ ÔÓÍÊÖÈÈ: [" << a->getXY()[0] << ", " << a->getXY()[1] << "] " << a->getBestScore() << "\n";
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
