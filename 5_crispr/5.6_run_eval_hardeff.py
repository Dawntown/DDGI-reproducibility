# %%
import os, sys
os.environ['CUDA_VISIBLE_DEVICES'] = '0'
import driver_genes as dg
dg.set_verbosity(False)
import pandas as pd
import torch

# %%
# num_deg_list = [10, 50, 100, 200, 500]
num_deg_list = sys.argv[1].split(',')
per_strength_list = [0.05, 0.25, 0.5, 0.75, 1]
extractor_type_list = ['attn']
# %%
for num_deg in num_deg_list:
    for per_strength in per_strength_list:
        for extractor_type in extractor_type_list:
            config_dict = dg.myutils.render_template(
                template_path='config_hardeff_template.yaml', 
                n_deg=num_deg, 
                per_strength=per_strength, 
                extractor_type=extractor_type
            )
            # train the model
            print(f'Training the {extractor_type} model with {num_deg} degrees and {per_strength} efficiency')
            pp = dg.Pipeline(config_dict=config_dict)
            pp.fit()
            pp.save_best_model()
            # reload the data with different efficiencies
            pp.args.dataset.adata_file = [
                '_'.join(pp.args.dataset.adata_file.split('_')[:-1]) + f'_{i}.h5ad'
                for i in per_strength_list
            ]
            pp.reload_data()
            # evaluate the model
            print(f'Evaluating the {extractor_type} model with {num_deg} degrees and all efficiencies')
            metrics_df, adata_pred = pp.evaluate(num_results=5, concat=False)
            metrics_df.to_csv(os.path.join(pp.args.trainer.output_dir, 'test_epoch_metrics_self.csv'), index=False)
            
            # evaluate with the true perturbation label
            metrics_func = dg.metrics.fetch_metrics(pp.args.trainer.test_metrics)
            os.makedirs(os.path.join(pp.args.trainer.output_dir, 'anndata'), exist_ok=True)
            metrics_list = []
            metrics_true_list = []
            for rep, ad_pred_dict in enumerate(adata_pred):
                for key, ad in ad_pred_dict.items():
                    ad_p = ad[ad.obs['perturbation'] != 'control']
                    cls = torch.from_numpy(ad_p.obsm['proba'].values)
                    psi = torch.from_numpy((ad_p.obs['perturbation'] == 'perturbed').values.astype(int).reshape(-1, 1))
                    psi_true = torch.from_numpy((
                        (ad_p.obs['perturbation'] == 'perturbed') &
                        (ad_p.obs['blurred'] == False)
                    ).values.astype(int).reshape(-1, 1))
                    # metrics for inpure label
                    metrics = dg.metrics.get_scalar_metrics_dict(metrics_func(cls, psi))
                    metrics['prefix'] = key
                    metrics['version'] = pp.args.trainer.version
                    metrics_list.append(metrics)
                    
                    # metrics for true label
                    metrics_true = dg.metrics.get_scalar_metrics_dict(metrics_func(cls, psi_true))
                    metrics_true['prefix'] = key
                    metrics_true['version'] = pp.args.trainer.version
                    metrics_true_list.append(metrics_true)
                    ad.write(os.path.join(
                        pp.args.trainer.output_dir, 
                        'anndata', 
                        f'{key}_rep{rep}.h5ad'
                    ))
            metrics_df = pd.DataFrame(metrics_list)
            metrics_df.to_csv(os.path.join(pp.args.trainer.output_dir, 'test_epoch_metrics.csv'), index=False)
            metrics_true_df = pd.DataFrame(metrics_true_list)
            metrics_true_df.to_csv(os.path.join(pp.args.trainer.output_dir, 'test_epoch_metrics_true.csv'), index=False)
                


# %%



