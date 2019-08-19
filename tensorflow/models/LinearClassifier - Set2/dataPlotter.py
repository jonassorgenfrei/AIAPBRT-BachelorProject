# Matplotlib
from matplotlib import pyplot as plt

def plot_features(validation_examples,validation_targets,training_examples,training_targets) :
    plt.figure(figsize=(13, 8))

    ax = plt.subplot(1, 2, 1)
    ax.set_title("Validation Data")

    ax.set_autoscaley_on(False)
    ax.set_ylim([32, 43])   #   limit y
    ax.set_autoscalex_on(False)
    ax.set_xlim([-126, -112])   # limit x
    plt.scatter(validation_examples["origin.x"],
                validation_examples["origin.y"],
                cmap="coolwarm",
                c=validation_targets["hit"])

    ax = plt.subplot(1,2,2)
    ax.set_title("Training Data")

    ax.set_autoscaley_on(False) 
    ax.set_ylim([32, 43])   # limit y
    ax.set_autoscalex_on(False)
    ax.set_xlim([-126, -112])   # limit x
    plt.scatter(training_examples["origin.x"],
                training_examples["origin.y"],
                cmap="coolwarm",
                c=training_targets["hit"])
    _ = plt.plot()
    plt.show()