import pickle
import prediction
import evaluate

with open('C:/Users/mz24b548/Desktop/scLinear_branch/gexp_fit_input','rb') as file:
        gex_matrix = pickle.load(file)

with open('C:/Users/mz24b548/Desktop/scLinear_branch/adt_fit_input','rb') as file:
        adt_matrix = pickle.load(file)

with open('C:/Users/mz24b548/Desktop/scLinear_branch/pipe','rb') as file:
        pipe = pickle.load(file)

with open('C:/Users/mz24b548/Desktop/scLinear_branch/gene_names','rb') as file:
        gene_names = pickle.load(file)

with open('C:/Users/mz24b548/Desktop/scLinear_branch/adt_names','rb') as file:
        adt_names = pickle.load(file)

with open('C:/Users/mz24b548/Desktop/scLinear_branch/evaluation_input_prediction','rb') as file:
        prediction = pickle.load(file)

with open('C:/Users/mz24b548/Desktop/scLinear_branch/evaluation_input_real','rb') as file:
        true_vals = pickle.load(file)

with open('C:/Users/mz24b548/Desktop/scLinear_branch/gexp_test_for_feature_importance','rb') as file:
        gexp_test = pickle.load(file)



pipe.fit(gex_matrix,adt_matrix,gex_names = gene_names,adt_names = adt_names,)
pipe.feature_importance(gexp_test)
# evaluate.evaluate(prediction,true_vals)
print('a')
