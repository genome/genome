-- Deploy web_task_params
-- requires: web_task

BEGIN;

CREATE TABLE IF NOT EXISTS web.task_params (
    genome_task_id character varying(255) NOT NULL,
    params text NOT NULL,
    CONSTRAINT task_params_pkey PRIMARY KEY (genome_task_id),
    CONSTRAINT task_params_genome_task_id_fkey FOREIGN KEY (genome_task_id) REFERENCES web.task(id)
);

COMMIT;
