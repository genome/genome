-- Deploy model_model_group_bridge
-- requires: model_model
-- requires: model_model_group

BEGIN;

CREATE TABLE IF NOT EXISTS model.model_group_bridge (
    model_group_id character varying(32) NOT NULL,
    model_id character varying(32) NOT NULL,
    CONSTRAINT model_group_bridge_pkey PRIMARY KEY (model_group_id, model_id),
    CONSTRAINT model_group_bridge_model_group_id_fkey FOREIGN KEY (model_group_id) REFERENCES model.model_group(id),
    CONSTRAINT model_group_bridge_model_id_fkey FOREIGN KEY (model_id) REFERENCES model.model(genome_model_id)
);

COMMIT;
