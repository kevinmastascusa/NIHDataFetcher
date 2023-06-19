"""
Disease Prevalence Analysis

Author: Kevin Mastascusa
Date: 6/19/2023

This program demonstrates data analysis and visualization of disease prevalence using NIH databases.
Provide the appropriate input to the functions for analysis and visualization.

Example Usage:
    analyze_disease_prevalence("asthma")

Output:
- Disease Prevalence Analysis:
    - The program fetches data from NIH databases related to the provided disease.
    - It performs analysis to determine the prevalence of the disease.
    - The output includes statistical information such as the number of cases, demographics, and prevalence rates.

TODO:
- Implement additional analysis functions based on specific disease-related queries.
- Customize visualizations and statistical analyses as needed.
"""

from Bio import Entrez


def analyze_disease_prevalence(disease_name: str):
    """
    Analyze disease prevalence using NIH databases.

    :param disease_name: The name of the disease.
    """
    Entrez.email = "kzm24@drexel.edu"  # Replace with your actual email

    # Perform data retrieval and analysis
    search_term = f"{disease_name}[MeSH Terms]"
    handle = Entrez.esearch(db="pubmed", term=search_term)
    record = Entrez.read(handle)
    handle.close()
    num_cases = record["Count"]

    # Example statistical analysis
    # You can customize this section based on your specific analysis requirements
    demographics = {"Age": ["0-10", "11-20", "21-30"], "Gender": ["Male", "Female", "Other"]}
    prevalence_rates = [0.2, 0.3, 0.1]

    # Output statistical information about disease prevalence
    print(f"Disease Prevalence Analysis for {disease_name}:")
    print(f"Number of Cases: {num_cases}")
    print("Demographics:")
    for demo, rates in zip(demographics.items(), prevalence_rates):
        print(f"- {demo[0]}: {', '.join(demo[1])}")
        print(f"  Prevalence Rate: {rates}")
    print("")


if __name__ == "__main__":
    # Example test call for the function
    disease_name = input("Enter the name of the disease: ")
    analyze_disease_prevalence(disease_name)
